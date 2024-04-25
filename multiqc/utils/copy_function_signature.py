"""
For strict type checking and auto-completion of functions like
multiqc.run(**kwargs), multiqc.parse_logs(**kwargs), multiqc.write_report(**kwargs)
that pass all keyword arguments to an inner function core.init_config(**kwargs).
Because those function are user-facing, we want the whole signature to be public.
Solution is borrowed entirely from the comment at
https://github.com/python/typing/issues/270#issuecomment-1346124813
"""

import sys

if sys.version_info < (3, 10):
    # ParamSpec and Concatenate are new in Python 3.10, so defining
    # no-op versions of decorators for earlier versions of Python.

    def copy_callable_signature(source):
        return lambda target: target

    def copy_method_signature(source):
        return lambda target: target

else:
    import functools
    from collections.abc import Callable
    from typing import Any, Concatenate, ParamSpec, TypeVar

    P = ParamSpec("P")
    T = TypeVar("T")

    def copy_callable_signature(source: Callable[P, T]) -> Callable[[Callable[..., T]], Callable[P, T]]:
        def wrapper(target: Callable[..., T]) -> Callable[P, T]:
            @functools.wraps(source)
            def wrapped(*args: P.args, **kwargs: P.kwargs) -> T:
                return target(*args, **kwargs)

            return wrapped

        return wrapper

    def copy_method_signature(
        source: Callable[Concatenate[Any, P], T],
    ) -> Callable[[Callable[..., T]], Callable[Concatenate[Any, P], T]]:
        def wrapper(target: Callable[..., T]) -> Callable[Concatenate[Any, P], T]:
            @functools.wraps(source)
            def wrapped(self: Any, /, *args: P.args, **kwargs: P.kwargs) -> T:
                return target(self, *args, **kwargs)

            return wrapped

        return wrapper
