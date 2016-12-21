#!/bin/bash

# Test runner for MultiQC.
if [ -z "${MULTIQC_TEST_ROOT:-}" ] ; then
    MULTIQC_TEST_ROOT="$(dirname $0)"/../../MultiQC_TestData
fi

#Canonicalize it
MULTIQC_TEST_ROOT="`readlink -f $MULTIQC_TEST_ROOT`"

if [ ! -d "$MULTIQC_TEST_ROOT/unit_tests" ] ; then
    cat <<.
Error: $MULTIQC_TEST_ROOT was not found or does not contain a
       unit_tests/ subdirectory.

To reduce bloat of the main GIT repository, the tests and test data are kept
in the MultiQC_TestData repository, which you need to check out from GitHub.

If the test data is somewhere other than the path above, set \$MULTIQC_TEST_ROOT
to the correct directory and retry.
.
    exit 1
fi

# You can run a specific test or all tests at once.
# Specific tests will run only in python3 by default, if you have it.


export RUN_SLOW_TESTS=${RUN_SLOW_TESTS:-0}

export PYTHONPATH="$(readlink -f "$(dirname $0)"/..)"
PY=python
if which python3 >/dev/null 2>&1 ; then
    PY=python3
fi

# All tests should run within the same CWD so they can easily find the files
# they need, without resorting to reading __FILE__.
cd "$MULTIQC_TEST_ROOT"

if [ "$*" == "" ] ; then
    python3 -munittest discover
    python  -munittest discover
else
    set -e
    for t in "$@" ; do $PY -munittest unit_tests.test_"$t" ; done
fi
