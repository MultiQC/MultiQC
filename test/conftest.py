# -*- coding: utf-8 -*-
# Copyright (C) 2017 by Per Unneberg
import pytest

def pytest_addoption(parser):
    group = parser.getgroup("multiqc", "multiqc test options")
    group.addoption("-M", "--multiqc-module", action="store",
                     default=False,
                     help="module to test",
                     dest="module")
