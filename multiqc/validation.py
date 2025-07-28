"""
Better validation of configs. Build on top of Pydantic, but prints more helpful error messages.
"""

import inspect
import logging
import re
from collections import defaultdict
from typing import Any, Dict, Set, Tuple

from PIL import ImageColor
from pydantic import BaseModel
from typeguard import TypeCheckError, check_type

from multiqc import config

logger = logging.getLogger(__name__)


class ModuleConfigValidationError(Exception):
    """
    Thrown in strict mode when some ValidatedConfig fails validation
    """

    def __init__(self, message: str, module_name: str):
        self.message = message
        self.module_name = module_name
        super().__init__()


_errors_by_cfg_path: Dict[Tuple[str, ...], Set[str]] = defaultdict(set)
_warnings_by_cfg_path: Dict[Tuple[str, ...], Set[str]] = defaultdict(set)


def get_current_module_name() -> str:
    callstack = inspect.stack()
    for n in callstack:
        if "multiqc/modules/" in n[1] and "base_module.py" not in n[1]:
            callpath = n[1].split("multiqc/modules/", 1)[-1]
            return f"{callpath}: "
    return ""


def add_validation_error(path_in_cfg: Tuple[str, ...], error: str):
    _errors_by_cfg_path[path_in_cfg].add(error)


def add_validation_warning(path_in_cfg: Tuple[str, ...], warning: str):
    _warnings_by_cfg_path[path_in_cfg].add(warning)


# To avoid cluttering stdout with duplicated messages
collapse_repeated_messages = False
_printed_errors: Set[str] = set()
_printed_warnings: Set[str] = set()


def reset():
    # Reset aggregating messages. Usually starts in the beginning of a module
    global _errors_by_cfg_path, _warnings_by_cfg_path, _printed_errors, _printed_warnings
    _errors_by_cfg_path = defaultdict(set)
    _warnings_by_cfg_path = defaultdict(set)
    _printed_errors = set()
    _printed_warnings = set()


def _print_warning(msg: str):
    if collapse_repeated_messages:
        if msg in _printed_warnings:
            return
        else:
            _printed_warnings.add(msg)
    logger.warning(msg)


def _print_error(msg: str):
    if collapse_repeated_messages:
        if msg in _printed_errors:
            return
        else:
            _printed_errors.add(msg)
    logger.error(msg)


def _print_large_dict(data: Dict[str, Any]) -> str:
    """
    Print a large dict - recursively find all large sub-dicts and sub-lists and print only first 10 items of each
    """

    def _truncate_recursive(obj: Any, depth: int = 0) -> str:
        if isinstance(obj, dict):
            items = list(obj.items())[:10]
            truncated = ", ".join(f"{k}: {_truncate_recursive(v, depth + 1)}" for k, v in items)
            if len(obj) > 10:
                truncated += f", ... ({len(obj) - 10} more items)"
            return f"{{{truncated}}}"
        elif isinstance(obj, (list, tuple, set)):
            items = list(obj)[:10]
            truncated = ", ".join(_truncate_recursive(x, depth + 1) for x in items)
            if len(obj) > 10:
                truncated += f", ... ({len(obj) - 10} more items)"
            return f"[{truncated}]"
        else:
            return str(obj)

    return _truncate_recursive(data)


class ValidatedConfig(BaseModel):
    """
    Wrapper of BaseModel with better validation error messages. Takes path_in_cfg to
    reference where exactly the error happened in error messages in nested configs
    """

    def __init__(self, path_in_cfg: Tuple[str, ...], **data: Any):
        # Format original data for concise error messages
        formatted_data = _print_large_dict(data)

        # Validate and cast values and remove invalid ones
        data = self.validate_fields(path_in_cfg, data)

        # Get all errors that start with path_in_cfg
        errs_by_cfg_path = {
            k: errs for k, errs in _errors_by_cfg_path.items() if k[: len(path_in_cfg)] == path_in_cfg and errs
        }
        wrns_by_cfg_path = {
            k: wrns for k, wrns in _warnings_by_cfg_path.items() if k[: len(path_in_cfg)] == path_in_cfg and wrns
        }
        if errs_by_cfg_path or wrns_by_cfg_path:
            modname = get_current_module_name()
            _path_str = ".".join(path_in_cfg)

            if wrns_by_cfg_path:
                msg = f"{modname}warnings while parsing {_path_str} {formatted_data}"
                if len(wrns_by_cfg_path) > 1:
                    msg += f" (total messages: {len(wrns_by_cfg_path)})"
                msg += ":"
                for path, _warnings in sorted(wrns_by_cfg_path.items()):
                    for _warning in _warnings:
                        msg += f"\n• '{path[-1]}': {_warning}"
                    _warnings_by_cfg_path[path].clear()  # Reset for interactive usage
                _print_warning(msg)

            if errs_by_cfg_path:
                msg = f"{modname}errors while parsing {_path_str} {formatted_data}"
                if len(errs_by_cfg_path) > 1:
                    msg += f" (total messages: {len(errs_by_cfg_path)})"
                msg += ":"
                for path, _errors in sorted(errs_by_cfg_path.items()):
                    for _error in _errors:
                        msg += f"\n• '{path[-1]}': {_error}"
                    _errors_by_cfg_path[path].clear()  # Reset for interactive usage
                _print_error(msg)

                if config.strict:
                    raise ModuleConfigValidationError(message=msg, module_name=modname)

        # By this point, data is a valid dict with only valid fields, but it still can
        # raise PydanticValidationError with unexpected errors
        super().__init__(**data)

    @classmethod
    def validate_fields(cls, path_in_cfg: Tuple[str, ...], values: Dict[str, Any]) -> Dict[str, Any]:
        # Remove underscores from field names (used for names matching reserved keywords, e.g. from_)
        for k in cls.model_fields.keys():
            if k.endswith("_") and k[:-1] in values:
                values[k] = values.pop(k[:-1])

        # Remove None values
        values = {k: v for k, v in values.items() if v is not None}

        # Check unrecognized fields
        filtered_values = {}
        available_fields = [k for k, v in cls.model_fields.items() if not v.deprecated]
        for name, val in values.items():
            if name not in cls.model_fields:
                add_validation_warning(
                    path_in_cfg + (name,),
                    f"unrecognized field. Available fields: {', '.join(available_fields)}",
                )
            else:
                filtered_values[name] = val
        values = filtered_values

        # Convert deprecated fields
        values_without_deprecateds: Dict[str, Any] = {}
        for name, val in values.items():
            new_name = cls.model_fields[name].deprecated
            if isinstance(new_name, str) and new_name not in values:
                add_validation_warning(path_in_cfg + (name,), f"deprecated field. Use '{new_name}' instead")
                values_without_deprecateds[new_name] = val
                continue
            values_without_deprecateds[name] = val
        values = values_without_deprecateds

        # Check missing fields
        for name, field in cls.model_fields.items():
            if field.is_required():
                if name not in values:
                    add_validation_error(path_in_cfg + (name,), "missing required field")
                    try:
                        values[name] = field.annotation() if field.annotation else None
                    except TypeError:
                        values[name] = None

        # Check types and validate specific fields
        corrected_values = {}
        for name, val in values.items():
            field = cls.model_fields[name]
            expected_type = field.annotation

            # Parse potential sub-models
            parse_method = getattr(cls, f"parse_{name}", None)
            if parse_method is not None:
                try:
                    val = parse_method(val, path_in_cfg=path_in_cfg + (name,))
                except Exception as e:
                    msg = f"failed to parse value '{val}': {e}"
                    add_validation_error(path_in_cfg + (name,), msg)
                    logger.debug(f"{msg}: {e}")
                    continue

            try:
                check_type(val, expected_type)
            except TypeCheckError as e:
                try:  # try casting to expected type?
                    if expected_type is not None:
                        if expected_type.__name__ in ["Optional", "Union"]:
                            expected_type = expected_type.__args__[0]
                    val = expected_type(val)  # type: ignore
                except Exception:
                    v_str = repr(val)
                    if len(v_str) > 20:
                        v_str = v_str[:20] + "..."
                    expected_type_str = str(expected_type).replace("typing.", "")
                    msg = rf"expected type '{expected_type_str}', got '{type(val).__name__}' {v_str}"
                    add_validation_error(path_in_cfg + (name,), msg)
                    logger.debug(f"{msg}: {e}")
                else:
                    corrected_values[name] = val
            else:
                corrected_values[name] = val

        values = corrected_values
        return values

    @classmethod
    def parse_color(
        cls,
        val,
        path_in_cfg: Tuple[str, ...],
    ):
        if val is None:
            return None
        if re.match(r"\d+,\s*\d+,\s*\d+", val):
            val_correct = f"rgb({val})"
        else:
            val_correct = val

        # Check if it's an rgba format - ImageColor.getrgb doesn't support rgba
        rgba_match = re.match(r"rgba\((\d+),\s*(\d+),\s*(\d+),\s*([0-9.]+)\)", val_correct)
        if rgba_match:
            # Validate that RGB values are in range 0-255
            try:
                r, g, b, a = rgba_match.groups()
                if not all(0 <= int(v) <= 255 for v in [r, g, b]):
                    raise ValueError("RGB values must be between 0 and 255")
                if not 0 <= float(a) <= 1:
                    raise ValueError("Alpha value must be between 0 and 1")
            except ValueError as e:
                add_validation_error(path_in_cfg, f"invalid color value '{val}': {e}")
                return None
            else:
                return val  # Return the original rgba string

        # For other formats, use ImageColor.getrgb
        try:
            ImageColor.getrgb(val_correct)
        except ValueError:
            add_validation_error(path_in_cfg, f"invalid color value '{val}'")
            return None
        else:
            return val
