from sequana.compare import RNADiffCompare
from sequana import sequana_data


def test_rnadiff_volcano():

    c = RNADiffCompare(
            sequana_data("rnadiff/rnadiff_onecond_1/tables/B3789-v1.surexpvsref.complete.xls"),
            sequana_data("rnadiff/rnadiff_onecond_2/tables/B0000-v1.surexpvsref.complete.xls"))


    c.plot_volcano()

    c.plot_common_major_counts("down")
    c.plot_common_major_counts("up")
    c.plot_common_major_counts("all")

    assert c.summary() == {"up1": 1295,   'up2': 56,
         'down1': 1325,
         'down2': 163,
         'common_down_r1_r2': 112,
         'common_up_r1_r2': 32}

