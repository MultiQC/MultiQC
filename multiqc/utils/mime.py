#!/usr/bin/env python

"""
File identification

Usage:
from multiqc.utils import mime
magic = mime.Magic()
with magic:
  modules = magic.guess_type(filepath)

"""

from tempfile import NamedTemporaryFile
import mimetypes
import subprocess
import os
import fnmatch
import re
import logging

from multiqc.utils import config


logger = logging.getLogger(__name__)


class MagicFile(object):
    """
    Identifies files using `file`
    """

    def __init__(self):
        super(MagicFile, self).__init__()

        self.magic_file = None
        self.magic_lines = []
        self.file_cmd = None

    def parse_contents(self, config, name):
        """
        Creates a `magic(5)` line from a (list of) exact patterns to match and
        adds it to the internal magic file.

        config: list of str or str
        name: handle to assign
        """
        if isinstance(config, str):
            config = [config]
        for line in config:
            flags = None
            if "    " in line:
                # reduce long spaces by half and tell magic to
                # require at least that much
                line = line.replace("  ", " ")
                flags = "W"
            line = line.replace(" ","\\ ")  # escape backslash
            line = line.replace("\t", "\\t")  # escape tab
            self.magic_lines.append(self.make_magic_line(pattern=line,
                                                         name=name, flags=flags))

    def make_magic_line(self, pattern, name,
                        offset=0, mode="search", len=4000, flags=None):
        """
        Creates a line in magic(5) format
        """
        return "{offset}\t{mode}/{len}{flags}\t{pattern}\t{name}\n".format(
            offset = offset,
            mode = mode,
            len = len,
            pattern = pattern,
            name = "multiqc/"+name,
            flags = "/"+flags if flags else ""
        )

    def __enter__(self):
        """
        Enter active mode with running `file` background process
        """
        logger.debug("Starting `file` background process")
        if self.file_cmd is None:
            self.magic_file = NamedTemporaryFile(mode="w")
            self.magic_file.writelines(self.magic_lines)
            self.magic_file.flush()
            self.file_cmd = self.start_file_cmd(self.magic_file.name)
        else:
            raise Exception("can't be in enter twice")
        return self

    def __exit__(self, a, b, c):
        """
        Exit active mode
        """
        logger.debug("Terminating `file` background process")
        self.file_cmd.stdin.close()
        self.file_cmd.wait()
        self.file_cmd = None
        logger.debug("Terminated `file` background process")

    def start_file_cmd(self, magicfile):
        """
        Start the `find` background process.

        magicfile: the magic(5) format file to use
        """
        cmd = [
            'file',
            '--magic-file', magicfile,
            '--no-buffer',  # no buffering so we can pipe lines
            '--brief',  # no repeating of filename in output
            '--files-from', '-', # read from stdin
        ]
        # export MAGIC to point to our magic file as well
        # this overrides the system magic file
        env = os.environ.copy()
        env['MAGIC'] = magicfile
        return subprocess.Popen(cmd, close_fds=True, shell=False,
                                stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                bufsize=1, universal_newlines=True, env=env)

    def guess_type_using_file(self, filepath):
        """
        Identifies the type of a file using the background `find` process.

        filepath: path to the file to identify
        returns: array of identifications (empty array on none)
        """
        if not self.file_cmd:
            # auto-wrap in `with` statement to run file cmd
            with self:
                return self.guess_type_using_file(filepath)

        logger.debug("  Matching file contents...")
        self.file_cmd.stdin.write(filepath+"\n")
        rval = [ item[8:]
                 for item in self.file_cmd.stdout.readline().split(",")
                 if item.startswith("multiqc/")]
        if len(rval) == 0:
            logger.debug("    No content matches found")
        else:
            logger.debug("    Content matched modules {}".format(rval))
        return rval


class MagicGlob(object):
    """
    Identifies files using globbing
    """
    def __init__(self):
        super(MagicGlob, self).__init__()
        self.globs = {}

    def parse_fn(self, config, name):
        """
        Adds a (list of) glob patterns to matching dictionary

        config: (list of) patterns
        name: module name
        """
        if isinstance(config, str):
            config = [config]
        for pattern in config:
            regex = re.compile(fnmatch.translate(pattern))
            if regex in self.globs:
                self.globs[regex] += [name]
            else:
                self.globs[regex] = [name]

    def guess_type_using_glob(self, filepath):
        """
        Matches a filepath to registered patterns
        (Only the first pattern is returned)

        filepath: path to check
        returns: array of matched modules
        """
        # FIXME: Sort regex dictionary, otherwise we may have unpredictable
        #        behavior depending on dict order.

        logger.debug("  Matching globs...")
        for regex in self.globs:
            if regex.match(filepath):
                logger.debug("    Glob matched {} from {}".format(repr(regex), self.globs[regex]))
                return self.globs[regex]
        logger.debug("    No globs matched")
        return []


class Magic(MagicFile, MagicGlob):
    def __init__(self):
        super(Magic, self).__init__()
        self.parse_config(config.sp)

    def parse_config(self, config, pre=""):
        for module in config:
            if module == 'fn':
                self.parse_fn(config[module], pre[:-1])
            elif module == 'contents':
                self.parse_contents(config[module], pre[:-1])
            else:
                self.parse_config(config[module], pre+module+"/")

    def guess_type(self, filepath):
        logger.debug("Checking \"{}\"".format(filepath))
        types = self.guess_type_using_file(filepath) + self.guess_type_using_glob(filepath)
        logger.debug("  Matched to modules {}".format(types))
        return types
