# -*- coding: utf-8 -*-
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2017 - Sequana Development Team
#
#  File author(s):
#      Thomas Cokelaer <thomas.cokelaer@pasteur.fr>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
import re

from sequana import logger
logger.name = __name__
import ruamel.yaml

__all__ = ["YamlDocParser"]


class YamlDocParser(object):
    """A simple parser to extract block content to be found in YAML files

    So as to create tooltips automatically in :ref:`sequanix`, one can comment
    YAML configuration file with block comments (see developers guide in
    :ref:`developers` )

    Once read and parsed, all block comments before top-level sections are to 
    be found in the dictionary :attr:`sections`.

    .. doctest::

        from sequana import snaketools
        from sequana.iotools import YamlDocParser
        module = snaketools.Module('quality_control')
        r = YamlDocParser(module.config)
        r.sections['fastqc']

    Those lines are removed from the docstring but available as a dictionary

    """
    def __init__(self, filename):
        """.. rubric:: constructor

        :param str filename: the YAML file to parse

        ::

            # main documentation

            # block comment
            section1:
                - item

            # block comment
            section2:

            # a comment

            section3:

        Here, section1 and section2 have block comments but not section3

        """
        self.filename = filename
        self.regex_section = re.compile("^[a-z,A-Z,_,0-9]+:")
        self._specials = ["choice__"]

        self.sections = {}
        self._read_data()
        self._parse_data()

    def _get_expected_sections(self):
        """Get the top level keys in the YAML file

        :return: list of top level sections' names"""
        with open(self.filename, "r") as fh:
            data = ruamel.yaml.load(fh.read(), ruamel.yaml.RoundTripLoader)
        keys = list(data.keys())
        return keys

    def _read_data(self):
        with open(self.filename, "r") as fh:
            self.data = fh.readlines()

    def _parse_data(self):
        """Parse the YAML file to get the block content (comments)
        before each top-level sections. See doc in the constructor

        Removes all # so that the block of comments can be interpreted as
        a standard docstring in Sequanix
        """
        current_block = []
        current_section = "docstring"

        # if we get a line that starts with #, this is a new comment or
        # part of a block comment. Otherwise, it means the current block
        # comment has ended.

        for this in self.data:
            # Beginning of a new section at top level
            if self.regex_section.findall(this):
                name = self.regex_section.findall(this)[0]
                current_section = name.strip(":")
                self.sections[current_section] = "".join(current_block)
                current_block = []
                current_section = None
            elif this.startswith('#'):    # a comment at top level
                current_block.append(this)
            elif this.strip() == "":      # an empty line
                #this was the main comment, or an isolated comment
                current_block = []
            else:  # a non-empty line to skip
                current_block = []

        for key in self._get_expected_sections():
            if key not in self.sections.keys():
                logger.warning("section %s not dealt by the parsing function" % key)

    def _block2docstring(self, section):
        if section not in self.sections.keys():
            logger.warning("%s not found in the yaml " % section)
            return
        comments = self.sections[section]
        docstring = []
        for line in comments.split("\n"):
            if "#############" in line:
                pass
            elif sum([this in line for this in self._specials]):
                pass
            else:
                if len(line)<2: # an empty line (to keep)
                    docstring.append("")
                else:
                    docstring.append(line[2:]) # strip the "# "characters
        docstring = "\n".join(docstring).strip()
        return docstring

    def _get_specials(self, section):
        """This method extracts data from the docstring

        Lines such as ::

            field_choice__ = ["a", "b"]

        are extracted. Where _choice is a special keyword to be
        found.

        """
        if section not in self.sections.keys():
            logger.warning("%s not found in the yaml " % section)
            return
        comments = self.sections[section]
        specials = {}
        for line in comments.split("\n"):
            if "#############" in line:
                pass
            elif sum([this in line for this in self._specials]):
                for special in self._specials:
                    line = line[2:]
                    key, value = line.split("=", 1)
                    key = key.strip().rstrip("__")
                    value = value.strip()
                    specials[key] = list(eval(value))
        return specials
