# Import -----------------------------------------------------------------------

import os
import pandas as pd
from reports import HTMLTable
from sequana.report_main import BaseReport
from sequana.report_submapping import SubMappingReport

# Class ------------------------------------------------------------------------
class MappingReport(BaseReport):
    """
    """
    def __init__(self, low_threshold=-3, high_threshold=3, 
            directory="report", **kargs):
        super(MappingReport, self).__init__(
                jinja_filename="mapping/index.html",
                directory=directory,
                output_filename="report_mapping.html", **kargs)
        self.jinja['title'] = "Mapping Report"
        self.low_t = low_threshold
        self.high_t = high_threshold

    def set_data(self, data):
        self.mapping = data

    def parse(self):
        self.mapping.plot_coverage(filename=self.directory + os.sep + 
                                            "coverage.png")
        
        self.mapping.plot_hist(filename=self.directory + os.sep + 
                                            "zscore_hist.png")

        df = pd.DataFrame()
        formatter = '<a target="_blank" alt={0} href="{1}">{0}</a>'
        for i in range(0, len(self.mapping), 500000):
            name = "mapping_{0}".format(i)
            stop = i + 500000
            if stop > len(self.mapping):
                stop = len(self.mapping)
            name = "{0}_{1}".format(name, stop)
            link = name + ".html"
            r = SubMappingReport(start=i, stop=stop, 
                    output_filename=name + ".html",
                    directory=self.directory)
            r.set_data(self.mapping)
            r.create_report()
            print(link)
            print(self.directory)
            df = df.append({"name": formatter.format(name, link)}, 
                ignore_index=True)
        html = HTMLTable(df)
        self.jinja['list_submapping'] = html.to_html(index=False)

        low_cov_df = self.mapping.get_low_coverage(self.low_t)
        merge_low_cov = self.mapping.merge_region(low_cov_df)
        html = HTMLTable(merge_low_cov)
        html.add_bgcolor("size")
        self.jinja['low_coverage'] = html.to_html(index=False)
        
        high_cov_df = self.mapping.get_high_coverage(self.high_t)
        merge_high_cov = self.mapping.merge_region(high_cov_df)
        html = HTMLTable(merge_high_cov)
        html.add_bgcolor("size")
        self.jinja['high_coverage'] = html.to_html(index=False)
