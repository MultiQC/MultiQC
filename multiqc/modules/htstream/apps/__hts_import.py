from inspect import isclass

# from pkgutil import iter_modules
from importlib import import_module

from os.path import dirname, basename, isfile, join
import glob


# global list of support
supported_apps = {}
modules = glob.glob(join(dirname(__file__), "*.py"))
__all__ = [basename(f)[:-3] for f in modules if isfile(f) and not f.startswith("__hts")]
mod_path = __name__[:-13]


# for (_, module_name, _) in iter_modules([package_dir]):
for module_name in __all__:
    if not module_name.startswith("__"):
        # import the module and assign it to a global variable
        module = import_module(f"{mod_path}.{module_name}")

        if module_name != "htstream_utils":
            for attribute_name in dir(module):
                attribute = getattr(module, attribute_name)

                if isclass(attribute):
                    # Add the class to this package's and reducer variables
                    supported_apps[attribute_name] = attribute

        else:
            supported_apps[module_name] = module
