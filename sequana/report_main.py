import os
import easydev
sepjoin = os.sep.join

from reports import Report
from sequana import version

import glob

def _get_template_path(name):
    # Is it a local directory ?
    if os.path.exists(name):
        return name
    else:
        main_path = easydev.get_package_location("sequana")
        return os.sep.join([main_path, 'sequana', "resources", "jinja", name])


class BaseReport(Report):
    """A Parent child for all reports created in Sequana


    """
    def __init__(self, jinja_filename, directory="report",     
                output_filename="test.html", **kargs):
        """.. rubric:: Constructor

        :param jinja_template: name of a directory (either local) or
            from sequana/share/templates where JINJA files are available. A file
            named index.html is required but may be renamed (with
            **output_filename** parameter).
        :param output_filename: name of the final HTML file.
        :param directory: name of the output directory (defaults to report)

        Parameters accepted by :class:`reports.Report` are also accepted.

        """
        # finds automatically the local directory or a directory to be found in
        # sequana distribution
        sequana_path = easydev.get_package_location('sequana')
        extra_css_path = sepjoin([sequana_path, "sequana", "resources", "css"])
        extra_css_list = glob.glob(extra_css_path + os.sep + "*css")

        searchpath = sepjoin([sequana_path, "sequana", "resources", "jinja"])

        super(BaseReport, self).__init__(searchpath, filename=output_filename, 
            template_filename=jinja_filename,
            directory=directory, extra_css_list=extra_css_list, **kargs)

        # This redefines the default name of the output (index.html) otherwise,
        # several reports will overwrite the default index.html. 
        self.filename = output_filename

        # That is where we will store all data to be used by the Jinja templates
        # in Report, everything is saved in jinja, but just to not get confused,
        # we will use data attribute for now
        #self.data = {}

        # Here, we defined default values from what is expected from the Jinja
        # template in share/templates/adapter_removal
        self.jinja['sequana_version'] = version

        # Common information to be filled (possibly)
        #self.data['command'] = "unset"
        self.jinja['dependencies'] =  self.get_table_dependencies('sequana').to_html()

    def parse(self):
        """populate the :attr:`data` attribute used by the JINJA templates

        """
        raise NotImplementedError

    def create_report(self, onweb=False):
        try:
            self.parse()
        except Exception as err:
            print(err)
            print("Parsing of %s failed" % self.input_filename)
        super(BaseReport, self).create_report(onweb=onweb)

    def read_snakefile(self, snakefile):
        """

        :param str filename:
        """
        with open(snakefile, "r") as fin:
            self.jinja['snakefile'] = fin.read()

    def read_configfile(self, configfile):
        with open(configfile, "r") as fin:
            self.jinja['config'] = fin.read()



class SequanaReport(BaseReport):
    def __init__(self, snakefile="Snakefile", configfile="config.yaml",
                 stats="stats.txt", **kargs):

        super(SequanaReport, self).__init__(
            jinja_filename="main/index.html", 
            directory="report",  
            output_filename="index.html",
            **kargs)

        try:
            self.read_snakefile(snakefile)
        except:
            pass

        try:
            self.read_configfile(configfile)
        except:
            pass

        try:
            from sequana.snakemake import SnakeMakeStats
            s = SnakeMakeStats(stats)
            s.plot()
            import pylab
            output_file = "snakemake_stats.png"
            pylab.savefig(self.directory + os.sep + output_file)
            self.jinja['snakemake_stats'] = output_file
        except:
            print('snakemake stats.txt not found')


        self.jinja['title'] = "Sequana Report"

    def parse(self):
        pass
 










