import tempfile
import os
from sequana.sphinxext import snakemakerule
from sequana.sphinxext import sequana_pipeline
from sphinx.application import Sphinx

data = """import sys, os
import sphinx
sys.path.insert(0, os.path.abspath('sphinxext'))
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    "sequana.sphinxext.snakemakerule",
    "sequana.sphinxext.sequana_pipeline"
    ]
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = "sequana"
copyright = "2016"
version = '1.0' 
release = "1.0"
exclude_patterns = []
add_module_names = False
pygments_style = 'sphinx'
intersphinx_mapping = {}
"""

def test_sequana_pipeline():
    res  = sequana_pipeline.get_rule_doc("quality_control")

    with tempfile.TemporaryDirectory() as tmpdir:

        # Create the conf and index in tmpdir
        with open(tmpdir+os.sep+"index.rst", "w") as fh:
            fh.write(".. snakemakerule:: dag\n")

        with open(tmpdir+os.sep+"conf.py", "w") as fh:
            print(fh.name)
            fh.write(data)

        app = Sphinx(tmpdir, tmpdir, tmpdir+"/temp", tmpdir, "html")
        app.build()

def test_doc():
    res  = snakemakerule.get_rule_doc("dag")
    res  = snakemakerule.get_rule_doc("fastqc_dynamic")

    try:
        res  = snakemakerule.get_rule_doc("dummy")
        assert False
    except FileNotFoundError:
        assert True
    except:
        assert False

    with tempfile.TemporaryDirectory() as tmpdir:

        # Create the conf and index in tmpdir
        with open(tmpdir+os.sep+"index.rst", "w") as fh:
            fh.write(".. snakemakerule:: dag\n")

        with open(tmpdir+os.sep+"conf.py", "w") as fh:
            print(fh.name)
            fh.write(data)

        # srcdir, confdir, outdir, doctreedir, buildername
        app = Sphinx(tmpdir, tmpdir, tmpdir+"/temp", tmpdir, "html")
        app.build()
