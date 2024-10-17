from sequana.lazy import pylab


def hist2D(
    x,
    y,
    data,
    hold=False,
    fontsize=12,
    grid=True,
    xlabel=None,
    ylabel=None,
    title=None,
    hist2d_dict={"bins": [100, 100], "contour": False, "norm": "log", "Nlevels": 6, "cmap": "BrBg"},
):

    from sequana.viz import Hist2D

    if hold is False:
        pylab.clf()

    data = data.loc[:, [x, y]].dropna()
    h = Hist2D(data)
    res = h.plot(**hist2d_dict)

    if xlabel:
        pylab.xlabel(xlabel, fontsize=fontsize)
    else:
        pylab.xlabel(x, fontsize=fontsize)

    if ylabel:
        pylab.ylabel(ylabel, fontsize=fontsize)
    else:
        pylab.ylabel(y, fontsize=fontsize)

    if title:
        pylab.title(title, fontsize=fontsize)

    if grid is True:
        pylab.grid(True)

    return pylab.gcf()
