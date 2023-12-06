import sequana.resources.canvas.bar as bar


def test_bar():

    data = [{"name": "A", "data": {"R1": 10, "R2": 90}}, {"name": "B", "data": {"R1": 90, "R2": 10}}]
    bar.stacked_bar("title", "ACGT", datalist=data)
