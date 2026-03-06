__all__ = [
    "ANOVA",
    "Boxplot",
    "Clustermap",
    "Corrplot",
    "Heatmap",
    "hinton",
    "Hist2D",
    "Imshow",
    "PCA",
    "BinaryPercentage",
    "ScatterHist",
    "plot_venn",
    "Volcano",
]

_module_map = {
    "ANOVA": "sequana.viz.anova",
    "Boxplot": "sequana.viz.boxplot",
    "Corrplot": "sequana.viz.corrplot",
    "Heatmap": "sequana.viz.heatmap",
    "Clustermap": "sequana.viz.heatmap",
    "hinton": "sequana.viz.hinton",
    "Hist2D": "sequana.viz.hist2d",
    "Imshow": "sequana.viz.imshow",
    "PCA": "sequana.viz.pca",
    "BinaryPercentage": "sequana.viz.plotly",
    "ScatterHist": "sequana.viz.scatter",
    "plot_venn": "sequana.viz.venn",
    "Volcano": "sequana.viz.volcano",
}


def __getattr__(name):
    if name in _module_map:
        import importlib

        module = importlib.import_module(_module_map[name])
        return getattr(module, name)
    raise AttributeError(f"module 'sequana.viz' has no attribute {name!r}")
