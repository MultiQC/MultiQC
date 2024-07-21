import inspect
import logging
import re
from collections import defaultdict
from typing import Dict, Set, Union, Optional, List

from pydantic import BaseModel, model_validator, ValidationError as PydanticValidationError
from typeguard import check_type, TypeCheckError
from PIL import ImageColor

from multiqc import config

logger = logging.getLogger(__name__)


class ConfigValidationError(Exception):
    def __init__(self, message: str, module_name: str):
        self.message = message
        self.module_name = module_name
        super().__init__()


_validation_errors_by_cls: Dict[str, Set[str]] = defaultdict(set)
_validation_warnings_by_cls: Dict[str, Set[str]] = defaultdict(set)


def add_validation_error(cls: Union[type, List[type]], error: str):
    cls_name = ".".join([c.__name__ for c in cls]) if isinstance(cls, list) else cls.__name__
    _validation_errors_by_cls[cls_name].add(error)


def add_validation_warning(cls: Union[type, List[type]], warning: str):
    cls_name = ".".join([c.__name__ for c in cls]) if isinstance(cls, list) else cls.__name__
    _validation_warnings_by_cls[cls_name].add(warning)


class ValidatedConfig(BaseModel):
    def __init__(self, _parent_class: Optional[type] = None, **data):
        _clss = [self.__class__]
        _cls_name = self.__class__.__name__
        _full_cls_name = _cls_name
        if _parent_class is not None:
            _clss = [_parent_class] + _clss
            _full_cls_name = f"{_parent_class.__name__}.{_cls_name}"

        try:
            super().__init__(**data, _clss=_clss)
        except PydanticValidationError:
            if not _validation_errors_by_cls.get(_cls_name) or not _validation_errors_by_cls.get(_full_cls_name):
                raise

        errors = _validation_errors_by_cls.get(_cls_name, set()) | _validation_errors_by_cls.get(_full_cls_name, set())
        warnings = _validation_warnings_by_cls.get(_cls_name, set()) | _validation_warnings_by_cls.get(
            _full_cls_name, set()
        )

        if errors or warnings:
            modname = ""
            callstack = inspect.stack()
            for n in callstack:
                if "multiqc/modules/" in n[1] and "base_module.py" not in n[1]:
                    callpath = n[1].split("multiqc/modules/", 1)[-1]
                    modname = f"{callpath}: "
                    break

            if warnings:
                logger.warning(f"{modname}Warnings while parsing {_full_cls_name} {data}:")
                for warning in sorted(warnings):
                    logger.warning(f"• {warning}")
                _validation_warnings_by_cls[_cls_name].clear()  # Reset for interactive usage
                _validation_warnings_by_cls[_full_cls_name].clear()  # Reset for interactive usage

            if errors:
                msg = f"{modname}Errors parsing {_full_cls_name} {data}"
                logger.error(msg)
                for error in sorted(errors):
                    logger.error(f"• {error}")
                    msg += f"\n• {error}"
                _validation_errors_by_cls[_cls_name].clear()  # Reset for interactive usage
                _validation_errors_by_cls[_full_cls_name].clear()  # Reset for interactive usage
                if config.strict:
                    raise ConfigValidationError(message=msg, module_name=modname)

    # noinspection PyNestedDecorators
    @model_validator(mode="before")
    @classmethod
    def validate_fields(cls, values):
        # Check unrecognized fields
        if not isinstance(values, dict):
            return values

        _clss = values.pop("_clss", [cls])

        # Remove underscores from field names (used for names matching reserved keywords, e.g. from_)
        for k, v in cls.model_fields.items():
            if k.endswith("_") and k[:-1] in values:
                values[k] = values.pop(k[:-1])

        # Check unrecognized fields
        filtered_values = {}
        for name, val in values.items():
            if name not in cls.model_fields:
                add_validation_warning(
                    _clss, f"unrecognized field '{name}'. Available fields: {', '.join(cls.model_fields.keys())}"
                )
            else:
                filtered_values[name] = val
        values = filtered_values

        # Convert deprecated fields
        values_without_deprecateds = {}
        for name, val in values.items():
            if cls.model_fields[name].deprecated:
                new_name = cls.model_fields[name].deprecated
                add_validation_warning(_clss, f"Deprecated field '{name}'. Use '{new_name}' instead")
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
                    val = parse_method(val)
                except Exception as e:
                    msg = f"'{name}': failed to parse value '{val}'"
                    add_validation_error(_clss, msg)
                    logger.debug(f"{msg}: {e}")
                    continue

            try:
                check_type(val, expected_type)
            except TypeCheckError as e:
                v_str = repr(val)
                if len(v_str) > 20:
                    v_str = v_str[:20] + "..."
                expected_type_str = str(expected_type).replace("typing.", "")
                msg = f"'{name}': expected type '{expected_type_str}', got '{type(val).__name__}' {v_str}"
                add_validation_error(_clss, msg)
                logger.debug(f"{msg}: {e}")
            else:
                corrected_values[name] = val

        values = corrected_values
        return values

    @classmethod
    def parse_color(cls, val):
        if val is None:
            return None
        if re.match(r"\d+,\s*\d+,\s*\d+", val):
            val_correct = f"rgb({val})"
        else:
            val_correct = val
        try:
            ImageColor.getrgb(val_correct)
        except ValueError:
            add_validation_error(cls, f"invalid color value '{val}'")
            return None
        else:
            return val
