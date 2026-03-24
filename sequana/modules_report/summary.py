#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Module to write summary.html have all information about the pipeline and
to visit other analysis"""
import importlib.metadata
import os

from sequana.lazy import pandas as pd
from sequana.modules_report.base_module import SequanaBaseModule
from sequana.utils import config
from sequana.utils.datatables_js import DataTable


class SummaryBase(SequanaBaseModule):
    def __init__(self, required_dir=None):
        super(SummaryBase, self).__init__(required_dir=required_dir)

    def dependencies(self):
        """Table with all python dependencies and a text file with tools
        needed and their versions.
        """
        html_table = self.get_table_dependencies()
        pypi = self.create_link("Pypi", "http://pypi.python.org")
        try:
            if "requirements" in self.json:
                req = self.copy_file(self.json["requirements"], "inputs")
            else:
                raise Exception
        except:
            try:
                req = self.json["requirements"]
            except:
                return
        req = self.create_link("requirements", req)
        content = "<p>Python dependencies (<b>{0}</b>){1}</p>".format(pypi, html_table, req)
        l, c = self.create_hide_section("Dep", "collapse/expand", content, hide=True)
        self.sections.append(
            {
                "name": "Dependencies {0}".format(self.add_float_right("<small>{0}</small>".format(l))),
                "anchor": "dependencies",
                "content": c,
            }
        )

    def get_table_dependencies(self):
        """Return dependencies of Sequana."""

        project_names = list()
        versions = list()
        links = list()
        pypi = "https://pypi.python.org/pypi/{0}"

        for dep in importlib.metadata.requires("sequana"):

            try:
                project_name, version = dep.split()
            except ValueError:
                project_name = dep
                version = "?"

            versions.append(version)
            project_names.append(project_name)
            links.append(pypi.format(project_name))

        df = pd.DataFrame({"package": project_names, "version": versions, "link": links})
        df["sort"] = df["package"].str.lower()
        df.sort_values(by="sort", axis=0, inplace=True)
        df.drop("sort", axis=1, inplace=True)
        datatable = DataTable(df, "dep")
        datatable.datatable.datatable_options = {
            "paging": "false",
            "bFilter": "false",
            "bInfo": "false",
            "bSort": "false",
        }
        datatable.datatable.set_links_to_column("link", "package")
        js = datatable.create_javascript_function()
        html = datatable.create_datatable()
        return js + "\n" + html


class SequanaReport(SummaryBase):
    """Write summary HTML report of an analysis. It contains all information
    about the pipeline used, input/output files and version of software.
    """

    def __init__(
        self,
        data,
        intro="",
        output_filename="summary.html",
        title="",
        workflow=True,
        Nsamples=1,
    ):
        """ """
        super(SequanaReport, self).__init__(required_dir=("js", "css"))
        self.json = data
        self.name = data.get("name", "undefined")
        self.title = f"Sequana Report Summary ({self.name})"
        self.wrappers = data.get("sequana_wrappers", "latest")
        self.intro = intro

        config.pipeline_version = data.get("pipeline_version", "latest")
        config.pipeline_name = data.get("name", "undefined")
        config.sequana_wrappers = self.wrappers

        self.create_report_content(workflow=workflow)
        self.create_html(output_filename)

    def create_report_content(self, workflow=True):
        """Create the report content."""
        self.sections = list()

        for section in config.summary_sections:
            self.sections.append(section)

        if workflow:
            self.workflow()
        self.dependencies()
        self.caller()

    def workflow(self):
        img = self.json["rulegraph"]
        dag_svg = self.include_svg_image(img, alt="workflow")
        snakefile = ".sequana/{}.rules".format(self.name)

        try:
            with open(snakefile, "r") as fp:
                code = self.add_code_section(fp.read(), "python")
            sf = self.create_hide_section("Sf", "Show/hide Snakemake file", code, hide=True)
            sf = "\n".join(sf)
        except IOError:
            sf = "no snakefile found in .sequana/"

        configfile = ".sequana/config.yaml"
        try:
            with open(configfile, "r") as fp:
                code = self.add_code_section(fp.read(), "yaml")
                c = self.create_hide_section("C", "Show/hide config file", code, hide=True)
                c = "\n".join(c)
        except IOError:
            c = "no config found in .sequana/"

        self.sections.append(
            {
                "name": "Workflow",
                "anchor": "workflow",
                "content": "<p>The following network shows the workflow of the pipeline. "
                "Blue boxes are clickable and redirect to dedicated reports."
                "</p>\n{0}\n"
                "<p>The analysis was performed with the following "
                '<a href="{3}">Snakemake</a> and <a href="{4}">configfile</a>:'
                "</p>\n"
                "<ul>\n"
                "    <li>{1}</li>\n"
                "    <li>{2}</li>\n"
                "</ul>".format(dag_svg, sf, c, snakefile, configfile),
            }
        )

    def get_table_versions(self):
        """Return third party tools from the requirements.txt and their versions."""

        # if no version.txt is found, return nothing
        if os.path.exists(".sequana/versions.txt") is False:
            return ""

        versions = []
        tools = []
        with open(".sequana/versions.txt", "r") as fin:
            for line in fin.readlines():
                # sometimes, if parsing if wrong, you may have more than 2 items...
                # e.g., warning in container that appear before expected output
                try:
                    tool, version = line.split()
                except Exception:
                    tool = line.split()[0]
                    version = "?"
                versions.append(version)
                tools.append(tool)

        df = pd.DataFrame({"tool": tools, "version": versions})
        try:
            df["sort"] = df["tool"].str.lower()
            df.sort_values(by="sort", axis=0, inplace=True)
            df.drop("sort", axis=1, inplace=True)
        except (KeyError, AttributeError):  # could be empty
            pass
        datatable = DataTable(df, "dep_and_version")
        datatable.datatable.datatable_options = {
            "paging": "false",
            "bFilter": "false",
            "bInfo": "false",
            "bSort": "false",
        }
        js = datatable.create_javascript_function()
        html = datatable.create_datatable()
        return js + "\n" + html

    def dependencies(self):
        """Table with all python dependencies and a text file with tools
        needed and their versions.
        """

        html_table_versions = self.get_table_versions()
        html_table_deps = self.get_table_dependencies()
        pypi = self.create_link("Pypi", "http://pypi.python.org")

        req = self.create_link("requirements", ".sequana/env.yml")

        content = (
            "<p>Third party tools can be found within containers (see config file abobe) if you use --use-apptainers option. Otherwise, here is a list of required dependencies and their versions.</p>"
            "<p>{3}</p>"
            "<p>Python dependencies (<b>{0}</b>){1}</p>".format(pypi, html_table_deps, req, html_table_versions)
        )
        l, c = self.create_hide_section("Dep", "collapse/expand", content, hide=True)
        self.sections.append(
            {
                "name": "Dependencies {0}".format(self.add_float_right("<small>{0}</small>".format(l))),
                "anchor": "dependencies",
                "content": c,
            }
        )

    def caller(self):
        try:
            with open(".sequana/info.txt") as fin:
                command = fin.readlines()
                command = "<pre>" + "\n".join([x for x in command if not x.startswith("#")]) + "</pre>"
        except Exception as err:
            print(err)
            command = "unknown"
        self.sections.append({"name": "Command", "anchor": "command", "content": command})
