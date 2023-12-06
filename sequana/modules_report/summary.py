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
import os

import easydev
from deprecated import deprecated

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
        content = (
            "<p>Dependencies downloaded from bioconda "
            "<b>{2}</b></p>"
            "<p>Python dependencies (<b>{0}</b>){1}</p>".format(pypi, html_table, req)
        )
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
        dep_list = easydev.get_dependencies("sequana")
        # if installed with conda, this will be empty
        if len(dep_list) == 0:
            return ""

        project_name = list()
        version = list()
        link = list()
        pypi = "https://pypi.python.org/pypi/{0}"
        for dep in dep_list:
            version.append(dep.version)
            project_name.append(dep.project_name)
            link.append(pypi.format(dep.project_name))
        df = pd.DataFrame({"package": project_name, "version": version, "link": link})
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


@deprecated(version="1.0", reason="will be removed in v1.0. Update your pipelines to use SequanaReport instead")
class SummaryModule(SummaryBase):
    """Write summary HTML report of an analysis. It contains all information
    about the pipeline used, input/output files and version of software.
    """

    def __init__(self, data, intro="", output_filename="summary.html"):
        """ """
        super().__init__()
        print("SummaryModule will be deprecated in the future. Use SummaryModule2 for now")
        self.json = data
        for this in ["inputs", "outputs"]:
            assert this in self.json
        self.title = "Sequana Report Summary"
        self.intro = intro
        self.create_report_content()
        self.create_html(output_filename)

    def create_report_content(self):
        """Create the report content."""
        self.sections = list()

        if self.json["inputs"]:
            self.pipeline_inputs()
        if self.json["outputs"]:
            self.pipeline_outputs()
        if self.json["html"]:
            self.pipeline_html()

        for section in config.summary_sections:
            self.sections.append(section)

        self.workflow()
        self.dependencies()

    def pipeline_inputs(self):
        """Links corresponding to the analysed input files."""
        # copy inputs in the input directory
        input_dir = "inputs"

        inputs = [self.copy_file(i, input_dir) for i in self.json["inputs"]]
        # create links list
        html_list = "<li>{0}</li>"
        links = [
            html_list.format(self.create_link(os.path.basename(i), i, newtab=False, download=True)) for i in inputs
        ]
        links = "<ul>{0}</ul>".format("\n".join(links))
        self.sections.append(
            {
                "name": "Inputs",
                "anchor": "input",
                "content": "<p>Link to the original data analysed.</p>\n" "{0}".format(links),
            }
        )

    def pipeline_outputs(self):
        """Links to important outputs generated by the pipeline"""
        # copy outputs in the output directory
        output_dir = "outputs"
        outputs = [self.copy_file(i, output_dir) for i in self.json["outputs"]]
        # create links list
        html_list = "<li>{0}</li>"
        links = [
            html_list.format(self.create_link(os.path.basename(i), i, newtab=False, download=True)) for i in outputs
        ]
        links = "<ul>{0}</ul>".format("\n".join(links))
        self.sections.append(
            {
                "name": "Outputs",
                "anchor": "outputs",
                "content": "<p>Link to the most important output files generated by the "
                "pipeline.</p>\n{0}".format(links),
            }
        )

    def pipeline_html(self):
        """Links to HTML pages created by the rules."""
        output_dir = "html"
        html = [self.copy_file(i, output_dir) for i in self.json["html"]]
        html_list = "<li>{0}</li>"
        links = [html_list.format(self.create_link(os.path.basename(i), i)) for i in html]
        links = "<ul>{0}</ul>".format("\n".join(links))
        self.sections.append(
            {
                "name": "External HTML",
                "anchor": "ext_html",
                "content": "<p>Link to HTML pages created by the pipeline.</p>\n{0}" "\n".format(links),
            }
        )

    def workflow(self):
        """Create the interactive DAG to navigate through pages."""
        snakefile = self.copy_file(self.json["snakefile"], "./inputs")
        configfile = self.copy_file(self.json["config"], "./inputs")
        # move the SVG file in the images directory
        img = self.copy_file(self.json["rulegraph"], "./images")
        dag_svg = self.include_svg_image(img, alt="rulegraph")
        with open(self.json["snakefile"], "r") as fp:
            code = self.add_code_section(fp.read(), "python")
        sf = self.create_hide_section("Sf", "Show/hide Snakemake file", code, hide=True)
        sf = "\n".join(sf)
        with open(self.json["config"], "r") as fp:
            code = self.add_code_section(fp.read(), "yaml")
        c = self.create_hide_section("C", "Show/hide config file", code, hide=True)
        c = "\n".join(c)
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


@deprecated(version="1.0", reason="will be removed in v1.0. Update your pipelines to use SequanaReport instead")
class SummaryModule2(SummaryBase):
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
        super(SummaryModule2, self).__init__(required_dir=("js", "css"))
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
                tool, version = line.split()
                versions.append(version)
                tools.append(tool)

        df = pd.DataFrame({"tool": tools, "version": versions})
        df["sort"] = df["tool"].str.lower()
        df.sort_values(by="sort", axis=0, inplace=True)
        df.drop("sort", axis=1, inplace=True)
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
            "<p>Third party tools can be found within containers (see config file) if you use --use-apptainers options. Otherwise, here is a list of required dependencies and their versions.</p>"
            "<p>{3}</p>"
            "<p>A conda environment was found and installed package are here: <b>{2}</b></p>"
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
                tool, version = line.split()
                versions.append(version)
                tools.append(tool)

        df = pd.DataFrame({"tool": tools, "version": versions})
        df["sort"] = df["tool"].str.lower()
        df.sort_values(by="sort", axis=0, inplace=True)
        df.drop("sort", axis=1, inplace=True)
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
            "<p>Third party tools can be found within containers (see config file) if you use --use-apptainers options. Otherwise, here is a list of required dependencies and their versions.</p>"
            "<p>{3}</p>"
            "<p>A conda environment was found and installed package are here: <b>{2}</b></p>"
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
