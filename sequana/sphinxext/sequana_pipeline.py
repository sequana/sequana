# coding: utf-8 -*-
"""snakemakerule

Defines a docutils directive for inserting simple docstring extracting from 
a snakemake rule (from sequana project).

::

    .. snakemakerule: dag

The name must be a valid sequana rule in the rules directory accesible via the
:class:`sequana.snaketools.Module` class

"""
import re
from docutils.parsers.rst import  directives
from docutils import nodes
from docutils.nodes import Body, Element
import urllib

def get_rule_doc(name):
    """Decode and return the docstring(s) of a sequana/snakemake rule."""

    url = "https://raw.githubusercontent.com/sequana/sequana_{}/master/README.rst".format(name)
    data = urllib.request.urlopen(url).read().decode("utf8")

    try:
        from sequana import Module
        m = Module(name)
        version = m.version
    except:
        version = "?"

    docstring = "**current version**:{}\n\n{}".format(version, data)

    return docstring


"""def get_rule_doc_snakemake(name):
    from sequana import Module
    import snakemake
    rule = Module(name)
    wf = snakemake.Workflow(rule.path + "/%s.rules" %  name)
    wf.include(rule.path + "/%s.rules" % name)
    docstring = list(wf.rules)[0].docstring
    return docstring
"""

class snakemake_base(Body, Element):
    def dont_traverse(self, *args, **kwargs):
        return []


class sequana_pipeline(snakemake_base):
    pass


def run(content, node_class, state, content_offset):
    node = node_class("")  # shall we add something here ?
    node.rule_docstring = get_rule_doc(content[0])
    state.nested_parse(content, content_offset, node)
    return [node]


def sequana_pipeline_directive(name, arguments, options, content, lineno,
                         content_offset, block_text, state, state_machine):
    return run(content, sequana_pipeline, state, content_offset)


def setup(app):
    #app.add_autodocumenter(RuleDocumenter)
    app.add_directive('sequana_pipeline', sequana_pipeline_directive, True, (0,0,0))

    # Add visit/depart methods to HTML-Translator:
    def visit_perform(self, node):
        # Ideally, we should use sphinx but this is a simple temporary solution
        from docutils import core
        from docutils.writers.html4css1 import Writer
        w = Writer()
        try:
            res = core.publish_parts(node.rule_docstring, writer=w)['html_body']
            self.body.append('<div class="">' + res + '</div>' )
            node.children = []
        except Exception as err:
            print(err)
            self.body.append('<div class=""> no docstring </div>' )


    def depart_perform(self, node):
        node.children = []

    def visit_ignore(self, node):
        node.children = []

    def depart_ignore(self, node):
        node.children = []

    app.add_node(sequana_pipeline,
                 html=(visit_perform, depart_perform),
                 latex=(visit_ignore, depart_ignore))







