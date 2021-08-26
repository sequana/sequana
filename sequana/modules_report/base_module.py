# coding: utf-8
#
#  This file is part of Sequana software
#
#  Copyright (c) 2016 - Sequana Development Team
#
#  File author(s):
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@pasteur.fr>,
#          <d.desvillechabrol@gmail.com>
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
"""Generic module is the parent module of all other module"""
import os
import shutil
import io
import base64

from jinja2 import Environment, PackageLoader

# from reports import HTMLTable
from sequana.utils import config
import colorlog

logger = colorlog.getLogger(__name__)


__all__ = ["SequanaBaseModule"]


class SequanaBaseModule(object):
    """Generic Module to write HTML reports.


    # to add a TOC, add this code::

        <div id="tocDiv">
        <ul id="tocList"> </ul>
        </div>


    """

    def __init__(self, template_fn="standard.html", required_dir=None):
        if required_dir is None:
            self.required_dir = ("css", "js", "images")
        else:
            self.required_dir = required_dir

        self.output_dir = config.output_dir
        self.path = "./"
        # Initiate jinja template
        env = Environment(loader=PackageLoader("sequana", "resources/templates/"))
        self.template = env.get_template(template_fn)
        self._init_report()
        self._fotorama_js_added = False

    def _init_report(self):
        """Create the report directory. All necessary directories are copied
        in working directory.
        """
        # Be aware of #465 issue. We need to check that the target file is
        # valid, in which case there is no need to copy the files.

        # Create report directory
        if os.path.isdir(config.output_dir) is False:
            os.mkdir(self.output_dir)

        for directory in self.required_dir:
            complete_directory = os.sep.join([self.output_dir, directory])
            if os.path.isdir(complete_directory) is False:
                os.mkdir(complete_directory)

        # Copy css/js necessary files
        for filename in config.css_list:
            target = os.sep.join([self.output_dir, "css"])
            if os.path.isfile(target) is False:
                shutil.copy(filename, target)

        for filename in config.js_list:
            target = os.sep.join([self.output_dir, "js"])
            if os.path.isfile(target) is False:
                shutil.copy(filename, target)

    def create_html(self, output_filename):
        """Create HTML file with Jinja2.

        :param str output_filename: HTML output filename
        """
        if output_filename is None:
            return
        report_output = self.template.render(config=config, module=self)
        with open(os.sep.join([config.output_dir, output_filename]), "w") as fp:
            print(report_output, file=fp)

    def create_link(self, name, target, newtab=True, download=False):
        """Create an HTML hyperlink with name and target.

        :param str target: the target url.
        :param bool newtab: open html page in a new tab.
        :param bool download: download the target.

        Return as string the HTML hyperlink to the target.
        """
        link = '<a href="{0}" '
        if newtab:
            link += 'target="_blank" '
        if download:
            link += 'download="{0}" '
        link += ">{1}</a>"
        return link.format(target, name)

    def create_hide_section(self, html_id, name, content, hide=False):
        """Create an hideable section.

        :param str html_id: add short id to connect all elements.
        :param str name: name of the hyperlink to hide or display the content.
        :param str content: hideable HTML content.
        :param bool hide: set if the first state is hiding or not.

        Return tuple that contains HTML hyperlink and hideable section.
        """
        link = "<a href='#1' class='show_hide{0}'>{1}</a>".format(html_id, name)
        content = "<div class='slidingDiv{0}'>\n{1}\n</div>".format(html_id, content)
        hidden = ""
        if hide:
            hidden = '\n$(".slidingDiv{0}").hide();'.format(html_id)
        js = """
<script type="text/javascript">
    $(document).ready(function(){{{1}
        $(".show_hide{0}").click(function(){{
            $(".slidingDiv{0}").slideToggle();
        }});
    }});
</script>
        """.format(
            html_id, hidden
        )
        content = js + content
        return link, content

    def copy_file(self, filename, target_dir):
        """Copy a file to a target directory in report dir. Return the
        relative path of your file.

        :param str filename: file to copy.
        :param str target_dir: directory where to copy.

        Return relative path of the new file location.
        """
        directory = config.output_dir + os.sep + target_dir
        try:
            os.makedirs(directory)
        except FileExistsError:
            if os.path.isdir(directory):
                pass
            else:
                msg = "{0} exist and it is not a directory".format(directory)
                logger.error(msg)
                raise FileExistsError
        try:
            shutil.copy(filename, directory)
        except FileNotFoundError:
            msg = "{0} doesn't exist".format(filename)
            raise FileNotFoundError(msg)
        return target_dir + os.sep + os.path.basename(filename)

    def add_float_right(self, content):
        """Align a content to right."""
        return '<div style="float:right">{0}</div>'.format(content)

    def add_code_section(self, content, language):
        """Add code in your html."""
        html = '<div class="code"><pre><code class="{0}">{1}' "</code></pre></div>"
        return html.format(language, content)

    def include_svg_image(self, filename, alt="undefined"):
        """Include SVG image in the html."""
        html = (
            '<object data="{0}" type="image/svg+xml">\n'
            '<img src="{0}" alt={1}></object>'
        )
        return html.format(filename, alt)

    def png_to_embedded_png(self, png, style=None, alt="", title=""):
        """Include a PNG file as embedded file."""
        import base64

        with open(png, "rb") as fp:
            png = base64.b64encode(fp.read()).decode()
        if style:
            html = '<img style="{0}" alt="{1}" title="{2}"'.format(style, alt, title)
        else:
            html = '<img alt="{}" title="{}"'.format(alt, title)
        return '{0} src="data:image/png;base64,{1}">'.format(html, png)

    def create_embedded_png(self, plot_function, input_arg, style=None, **kwargs):
        """Take as a plot function as input and create a html embedded png
        image. You must set the arguments name for the output to connect
        buffer.
        """
        buf = io.BytesIO()
        # add buffer as output of the plot function
        kwargs = dict({input_arg: buf}, **kwargs)
        try:
            plot_function(**kwargs)
            html = "<img "
            if style:
                html += 'style="{0}"'.format(style)
            html += 'src="data:image/png;base64,{0}"/>'.format(
                base64.b64encode(buf.getvalue()).decode("utf-8")
            )
            buf.close()
        except Exception as err:
            print(err)
            html = "image not created"
        return html

    def create_combobox(self, path_list, html_id, newtab=True):
        """Create a dropdown menu with QueryJS.

        :param list path_list: list of links.

        return html div and js script as string.
        """
        option_list = (
            "<li>{0}</li>\n".format(
                self.create_link(os.path.basename(path), path, newtab)
            )
            for path in path_list
        )
        html = """
<div id="jq-dropdown-{1}" class="jq-dropdown jq-dropdown-tip jq-dropdown-scroll">
    <ul class="jq-dropdown-menu">
        {0}
    </ul>
</div>
<a href="#" data-jq-dropdown="#jq-dropdown-{1}">Subchromosome</a>
        """.format(
            "\n".join(option_list), html_id
        )
        return html

    def add_fotorama(
        self,
        files,
        width=600,
        height=800,
        loop=True,
        thumbnails=True,
        file_thumbnails=None,
        captions=None,
    ):

        if self._fotorama_js_added is False:
            script = """
        <!-- jQuery 1.8 or later, 33 KB -->
        <!--<script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.1/jquery.min.js"></script>-->
        <!-- Fotorama from CDNJS, 19 KB -->
        <link  href="https://cdnjs.cloudflare.com/ajax/libs/fotorama/4.6.4/fotorama.css" rel="stylesheet">
        <script src="https://cdnjs.cloudflare.com/ajax/libs/fotorama/4.6.4/fotorama.js"></script>
        """
            self._fotorama_js_added = True
        else:
            script = ""

        if captions:
            if len(files) != len(captions):
                raise ValueError(
                    "captions and files must be of same length with 1-to-1 mapping"
                )
        else:
            captions = [filename.split("/")[-1] for filename in files]

        script += '<div class="fotorama" fzyz-keyboard="true" '
        if thumbnails is True:
            script += ' data-nav="thumbs"'
        if loop is True:
            script += ' data-loop="true"'
        script += ' data-width="{}"  data-height="{}"'.format(width, height)
        script += ">"
        for filename, caption in zip(files, captions):
            script += '<img src="{}" data-caption="{}">'.format(filename, caption)
        script += "</div>"
        return script
