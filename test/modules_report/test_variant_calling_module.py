from sequana.modules_report.variant_calling import VariantCallingModule
from sequana.utils import config

from .. import test_dir


def test_variant_calling_module(tmpdir):
    directory = tmpdir.mkdir("test_variant_calling_module")
    config.output_dir = str(directory)
    config.sample_name = "JB409847"
    VariantCallingModule(f"{test_dir}/data/csv/JB409847.vc.csv")
