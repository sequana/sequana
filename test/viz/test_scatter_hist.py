from sequana.viz import ScatterHist


def test1():

    s = ScatterHist(x=[1, 2, 3, 4], y=[3, 4, 5, 6])
    s.plot()
    s.plot(scatter_position="top right")
    s.plot(scatter_position="top left")
    s.plot(scatter_position="bottom right")
