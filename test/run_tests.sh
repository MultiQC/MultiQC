#!/bin/bash

# Simple test wrapper.
# You can run a specific test or all tests at once.
# Specific tests will run only in python3 by default, if you have it.

cd "`dirname $0`"/..

export RUN_SLOW_TESTS=${RUN_SLOW_TESTS:-0}

PY=python
if which python3 >/dev/null 2>&1 ; then
    PY=python3
fi

if [ "$*" == "" ] ; then
    python3 -munittest discover
    python  -munittest discover
else
    set -e
    for t in "$@" ; do $PY -munittest test.test_"$t" ; done
fi
