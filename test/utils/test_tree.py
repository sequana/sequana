from sequana.utils import tree


def test_tree():
    t = tree.HTMLDirectory()
    t.get_html()
    t = tree.HTMLDirectory(pattern=".py")
    t.get_html()
