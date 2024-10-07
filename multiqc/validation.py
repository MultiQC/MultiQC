import inspect
import logging
import re
from collections import defaultdict
from typing import Any, Dict, List, Optional, Set, Union, Type, cast

from PIL import ImageColor
from pydantic import BaseModel, model_validator
from pydantic import ValidationError as PydanticValidationError
from typeguard import TypeCheckError, check_type

from multiqc import config

logger = logging.getLogger(__name__)


class ModuleConfigValidationError(Exception):
    # Thrown in strict mode when some ValidatedConfig fails validation
    def __init__(self, message: str, module_name: str):
        self.message = message
        self.module_name = module_name
        super().__init__()


_errors_by_cls: Dict[str, Set[str]] = defaultdict(set)
_warnings_by_cls: Dict[str, Set[str]] = defaultdict(set)


def get_current_module_name() -> str:
    callstack = inspect.stack()
    for n in callstack:
        if "multiqc/modules/" in n[1] and "base_module.py" not in n[1]:
            callpath = n[1].split("multiqc/modules/", 1)[-1]
            return f"{callpath}: "
    return ""


def add_validation_error(cls: Union[type, List[type]], error: str):
    cls_name = ".".join([c.__name__ for c in cls]) if isinstance(cls, list) else cls.__name__
    _errors_by_cls[cls_name].add(error)


def add_validation_warning(cls: Union[type, List[type]], warning: str):
    cls_name = ".".join([c.__name__ for c in cls]) if isinstance(cls, list) else cls.__name__
    _warnings_by_cls[cls_name].add(warning)


# To avoid cluttering stdout with duplicated messages
collapse_repeated_messages = False
_printed_errors: Set[str] = set()
_printed_warnings: Set[str] = set()


def reset():
    # Reset aggregating messages. Usually starts in the beginning of a module
    global _errors_by_cls, _warnings_by_cls, _printed_errors, _printed_warnings
    _errors_by_cls = defaultdict(set)
    _warnings_by_cls = defaultdict(set)
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


class ValidatedConfig(BaseModel):
    def __init__(self, _clss: Optional[List[Type["ValidatedConfig"]]] = None, **data: Any):
        _cls_name = self.__class__.__name__
        _classes: List[Type["ValidatedConfig"]] = []
        if _clss:
            for _c in _clss:
                _classes.append(_c)
        _classes.append(self.__class__)
        _full_cls_name = ".".join(_c.__name__ for _c in _classes)

        try:
            super().__init__(**data, _clss=_clss)
        except PydanticValidationError:
            if not _errors_by_cls.get(_cls_name) or not _errors_by_cls.get(_full_cls_name):
                # Unhandled PydanticValidationError
                raise

        errors = _errors_by_cls.get(_cls_name, set()) | _errors_by_cls.get(_full_cls_name, set())
        warnings = _warnings_by_cls.get(_cls_name, set()) | _warnings_by_cls.get(_full_cls_name, set())
        if errors or warnings:
            modname = get_current_module_name()

            if warnings:
                msg = f"{modname}warnings while parsing {_full_cls_name} {data} (total: {len(warnings)}):"
                for warning in sorted(warnings):
                    msg += f"\n• {warning}"
                _print_warning(msg)

                # Reset for interactive usage
                _warnings_by_cls[_cls_name].clear()
                _warnings_by_cls[_full_cls_name].clear()

            if errors:
                msg = f"{modname}errors while parsing {_full_cls_name} {data} (total: {len(errors)}):"
                for err in sorted(errors):
                    msg += f"\n• {err}"
                _print_error(msg)

                # Reset for interactive usage
                _errors_by_cls[_cls_name].clear()
                _errors_by_cls[_full_cls_name].clear()

                if config.strict:
                    raise ModuleConfigValidationError(message=msg, module_name=modname)

    # noinspection PyNestedDecorators
    @model_validator(mode="before")
    @classmethod
    def validate_fields(cls, values: Any) -> Dict[str, Any]:
        # Check unrecognized fields
        if not isinstance(values, dict):
            return values
        values = cast(Dict[str, Any], values)
        _clss = cast(List[Type[ValidatedConfig]], values.pop("_clss", None) or [cls])

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
                    _clss,
                    f"unrecognized field '{name}'. Available fields: {', '.join(available_fields)}",
                )
            else:
                filtered_values[name] = val
        values = filtered_values

        # Convert deprecated fields
        values_without_deprecateds = {}
        for name, val in values.items():
            if cls.model_fields[name].deprecated:
                new_name = cls.model_fields[name].deprecated
                add_validation_warning(_clss, f"deprecated field '{name}'. Use '{new_name}' instead")
                if new_name not in values:
                    values_without_deprecateds[new_name] = val
            else:
                values_without_deprecateds[name] = val
        values = values_without_deprecateds

        # Check missing fields
        for name, field in cls.model_fields.items():
            if field.is_required():
                if name not in values:
                    add_validation_error(_clss, f"missing required field '{name}'")
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
                    val = parse_method(val, _clss=_clss)
                except Exception as e:
                    msg = f"'{name}': failed to parse value '{val}'"
                    add_validation_error(_clss, msg)
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
                    msg = f"'{name}': expected type '{expected_type_str}', got '{type(val).__name__}' {v_str}"
                    add_validation_error(_clss, msg)
                    logger.debug(f"{msg}: {e}")
                else:
                    corrected_values[name] = val
            else:
                corrected_values[name] = val

        values = corrected_values
        return values

    @classmethod
    def parse_color(cls, val, _clss=None):
        if val is None:
            return None
        if re.match(r"\d+,\s*\d+,\s*\d+", val):
            val_correct = f"rgb({val})"
        else:
            val_correct = val
        try:
            ImageColor.getrgb(val_correct)
        except ValueError:
            add_validation_error(_clss or cls, f"invalid color value '{val}'")
            return None
        else:
            return val
