from sequana import sequana_config_path, sequana_data


def test_config_directory():
    assert sequana_config_path


def test_sequana_data():

    try:
        sequana_data()
        assert False
    except ValueError:
        assert True
    except:
        assert False

    sequana_data("Hm2_GTGAAA_L005_R1_001.fastq.gz", "data")
