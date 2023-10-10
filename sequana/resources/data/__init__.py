"""

Some useful data sets to be used in the analysis


The command :func:`sequana.sequana_data` may be used to retrieved data from 
this package. For example, a small but standard reference (phix) is used in 
some NGS experiments. The file is small enough that it is provided within
sequana and its filename (full path) can be retrieved as follows::

    from sequana import sequana_data
    fullpath = sequana_data("phiX174.fa", "data")

Other files stored in this directory will be documented here.



"""


#: List of adapters used in various sequencing platforms
adapters = {
    "adapters_netflex_pcr_free_1_fwd": "adapters_netflex_pcr_free_1_fwd.fa",
    "adapters_netflex_pcr_free_1_rev": "adapters_netflex_pcr_free_1_rev.fa",
}
