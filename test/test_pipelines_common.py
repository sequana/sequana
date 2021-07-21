from easydev import AttrDict
import argparse

from sequana import sequana_data





def test_print_version():
    from sequana.pipelines_common import print_version
    try:print_version("quality_control")
    except:pass
    try:print_version("sequana_dummy")
    except:pass

def test_get_pipeline_location():
    from sequana.pipelines_common import get_pipeline_location
    get_pipeline_location("quality_control")


