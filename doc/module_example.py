import pandas as pd
from numpy import random

from sequana.modules_report.base_module import SequanaBaseModule
from sequana.utils.datatables_js import DataTable


class MyModule(SequanaBaseModule):
    def __init__(self, df, output="mytest.html"):
        super().__init__()
        self.data = df
        self.summary = self.data.describe().to_frame()

        self.title = "Super Module"
        self.create_report_content()
        self.create_html(output)

    def create_report_content(self):
        self.sections = list()
        self.add_table()
        self.add_image()

    def add_table(self):
        df = self.summary.copy()
        df.columns = ["data"]
        df["url"] = ["http://sequana.readthedocs.org"] * len(df)

        table = DataTable(df, "table", index=True)
        table.datatable.datatable_options = {
            "scrollX": "300px",
            "pageLength": 15,
            "scrollCollapse": "true",
            "dom": "tB",
            "paging": "false",
            "buttons": ["copy", "csv"],
        }
        table.datatable.set_links_to_column("url", "data")

        js = table.create_javascript_function()
        html_tab = table.create_datatable(float_format="%.3g")
        html = "{} {}".format(html_tab, js)

        self.sections.append({"name": "Table", "anchor": "table", "content": html})

    def add_image(self):
        import pylab

        def plotter(filename):
            pylab.ioff()
            self.data.hist()
            pylab.savefig(filename)

        html = self.create_embedded_png(plotter, "filename", style="width:65%")
        self.sections.append({"name": "Image", "anchor": "table", "content": html})


# Let us create some data.
df = pd.Series(random.randn(10000))

# and pass it as a first argument.
MyModule(df, "report_example.html")
