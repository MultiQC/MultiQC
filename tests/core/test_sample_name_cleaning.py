#!/usr/bin/python

"""Test that the sample cleaning logic works as expected."""

import unittest
from copy import copy

from multiqc import config
from multiqc.base_module import BaseMultiqcModule

# Capture default fn_clean_exts, fn_clean_trim before we trash them.
def_fn_clean_exts = copy(config.fn_clean_exts)
def_fn_clean_trim = copy(config.fn_clean_trim)


class T(unittest.TestCase):
    def setUp(self):
        self.bmm = BaseMultiqcModule()
        self.scrub_defaults = False
        self.config = {}

    def test_noop(self):
        self.scrub_defaults = True

        self.assertClean("foo.bar.fastq.gz", "foo.bar.fastq.gz")

    def test_default_trim(self):
        self.assertClean("foo.bar.fastq.gz", "foo.bar")

    def test_custom_clean_ext(self):
        self.config = dict(extra_fn_clean_exts=["chop_this_off"])

        self.assertClean("foo.chop_this_off", "foo")

    def test_regex_keep(self):
        self.config = dict(extra_fn_clean_exts=[{"type": "regex_keep", "pattern": "abc..X"}])

        self.assertClean("foo_abc12X_bar", "abc12X")

        self.assertClean("foo_abc123_bar", "foo_abc123_bar")

    # Helper functions

    def assertClean(self, unclean, clean):
        self.assertEqual(self.clean(unclean), clean)

    def clean(self, unclean, root=None):
        """Refresh the configuration, then invoke the cleanup logic with
        self.config upon the input sample name.
        """
        if self.scrub_defaults:
            config.fn_clean_exts[:] = []
            config.fn_clean_trim[:] = []
        else:
            config.fn_clean_exts[:] = def_fn_clean_exts
            config.fn_clean_trim[:] = def_fn_clean_trim
        config._add_config(self.config, None)

        return self.bmm.clean_s_name(unclean, root)


if __name__ == "__main__":
    unittest.main()
