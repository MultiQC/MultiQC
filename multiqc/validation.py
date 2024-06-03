import inspect
import logging

from pydantic import BaseModel, ValidationError, model_validator
from typeguard import check_type, TypeCheckError


logger = logging.getLogger(__name__)


class ConfigValidationError(Exception):
    def __init__(self, module_name: str):
        self.module_name = module_name
        super().__init__()


_validation_errors = []


class ValidatedConfig(BaseModel):
    def __init__(self, **data):
        try:
            super().__init__(**data)
        except ValidationError as e:
            if not _validation_errors:
                raise
            else:
                # errors are already added into plot.pconfig_validation_errors by a custom validator
                logger.debug(e)

        if _validation_errors:
            # Get module name
            modname = ""
            callstack = inspect.stack()
            for n in callstack:
                if "multiqc/modules/" in n[1] and "base_module.py" not in n[1]:
                    callpath = n[1].split("multiqc/modules/", 1)[-1]
                    modname = f"{callpath}: "
                    break

            plot_type = self.__class__.__name__.replace("Config", "")
            logger.error(f"{modname}Invalid {plot_type} plot configuration {data}:")
            for error in _validation_errors:
                logger.error(f"• {error}")
            _validation_errors.clear()  # Reset for interactive usage
            raise ConfigValidationError(module_name=modname)

    # noinspection PyNestedDecorators
    @model_validator(mode="before")
    @classmethod
    def validate_fields(cls, values):
        # Check unrecognized fields
        filtered_values = {}
        for name, val in values.items():
            if name not in cls.model_fields:
                _validation_errors.append(
                    f"unrecognized field '{name}'. Available fields: {', '.join(cls.model_fields.keys())}"
                )
            else:
                filtered_values[name] = val
        values = filtered_values

        # Convert deprecated fields
        values_without_deprecateds = {}
        for name, val in values.items():
            if cls.model_fields[name].deprecated:
                new_name = cls.model_fields[name].deprecated
                logger.debug(f"Deprecated field '{name}'. Use '{new_name}' instead")
                if new_name not in values:
                    values_without_deprecateds[new_name] = val
            else:
                values_without_deprecateds[name] = val
        values = values_without_deprecateds

        # Check missing fields
        for name, field in cls.model_fields.items():
            if field.is_required():
                if name not in values:
                    _validation_errors.append(f"missing required field '{name}'")

        # Check types
        for name, val in values.items():
            field = cls.model_fields[name]
            expected_type = field.annotation
            try:
                check_type(val, expected_type)
            except TypeCheckError as e:
                v_str = repr(val)
                if len(v_str) > 20:
                    v_str = v_str[:20] + "..."
                expected_type_str = str(expected_type).replace("typing.", "")
                msg = f"'{name}': expected type '{expected_type_str}', got '{type(val).__name__}' {v_str}"
                _validation_errors.append(msg)
                logger.debug(f"{msg}: {e}")

        return values
