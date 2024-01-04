def stacked_bar(tag, title, datalist):
    """

    data list should be a list of dictionary formatted as follows

    {"name": "A"
     "data": {
        "R1_mapped": 50,
        "R2_mapped": 50,
        "R1_unmapped": 50,
        "R2_unmapped": 50,
        }
    }
    """
    dataitems = ""
    for item in datalist:
        datatext = []
        for k, v in item["data"].items():
            datatext.append('{y:%s,label:"%s"}' % (v, k))
        datatext = ",\n        ".join(datatext)

        params = {"name": item["name"], "datatext": datatext}

        dataitems += (
            """
      {
        type: "stackedBar100",
        showInLegend: true,
        name: "%(name)s",
        dataPoints: [
        %(datatext)s
        ]
        },
    """
            % params
        )

    metadata = {"tag": tag, "title": title, "dataitems": dataitems}

    script = """
<script type="text/javascript">
  window.onload = function () {
    var chart = new CanvasJS.Chart("chartContainer%(tag)s",
    {
      theme: "theme2",
      title:{
        text: "%(title)s"
      },
      animationEnabled: true,
      axisY:{
        title: "percent"
      },
      legend :{
        horizontalAlign: 'center',
        verticalAlign: 'bottom'
      },
      toolTip: {
        shared: true
      },
      data:[
        %(dataitems)s
      ]

    });

chart.render();
}
</script>
"""
    return script % metadata


class CanvasBar(object):
    """ """

    def __init__(self, data, title="", tag="", xlabel="", ylabel_max_length=20, links=None, **kargs):
        """

        :param data: a dataframe with name, value, url
        """
        self.title = title
        self.tag = tag

        dataitems = "["
        for x, y in data.iterrows():
            formatter = '{y:%(value)s,label:"%(name)s",click:onClick,url:"%(url)s"},'
            dataitems += formatter % y
        dataitems += "]"

        self.metadata = {"dataitems": dataitems, "title": self.title, "xlabel": xlabel, "tag": self.tag}

    def to_html(self, options={"maxrange": None}):
        script = """
    var chart = new CanvasJS.Chart("chartContainer%(tag)s",
    {
      theme: "theme2",
      title:{
        text: "%(title)s"
      },
      animationEnabled: true,
      axisY:{
        title: "%(xlabel)s",
        %(axisYmaximum)s
      },
      exportFileName: "summary",
      exportEnabled: true,
      legend :{
        horizontalAlign: 'center',
        verticalAlign: 'bottom'
      },
      toolTip: {
        shared: true
      },
      data:[
      {
          type: "bar",
          dataPoints: %(dataitems)s
      },
      ]
    });

chart.render();

"""

        params = self.metadata.copy()
        if options["maxrange"]:
            params["axisYmaximum"] = "maximum:%s," % options["maxrange"]
        else:
            params["axisYmaximum"] = ""
        return script % params


# <div id="chartContainer%(tag)s" style="height: 300px; width: 100%%;">
