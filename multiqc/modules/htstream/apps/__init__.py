from inspect import isclass

# from pkgutil import iter_modules
from pathlib import Path
from importlib import import_module

from os.path import dirname, basename, isfile, join
import glob

# global list of support
globals()["supported_apps"] = []

modules = glob.glob(join(dirname(__file__), "*.py"))
__all__ = [basename(f)[:-3] for f in modules if isfile(f) and not f.endswith("__init__.py")]


# for (_, module_name, _) in iter_modules([package_dir]):
for module_name in __all__:

    # import the module and assign it to a global variable
    module = import_module(f"{__name__}.{module_name}")

    # add module name to list of support modules
    globals()["supported_apps"].append(module_name)

    if module_name != "htstream_utils":

        for attribute_name in dir(module):

            attribute = getattr(module, attribute_name)

            if isclass(attribute):

                # Add the class to this package's and reducer variables
                globals()[attribute_name] = attribute

    else:
        globals()[module_name] = module
