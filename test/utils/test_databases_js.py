



import pandas as pd
from sequana.utils.datatables_js import DataTable


def test_df():
    df = pd.DataFrame([[1,2], [3,4], ['a','b']])
    df = df.T
    df.columns = ['a','b','c']


    datatable = DataTable(df, 'data')


    # rerun constructor to have the tooltips
    datatable = DataTable(df, 'data')
    datatable.datatable_options = {'pageLength': 15,
                                   'dom': 'Bfrtip',
                                   'buttons': ['copy', 'csv']}

    datatable.datatable.set_tooltips_to_column('c','a')
    js = datatable.create_javascript_function()
    html = datatable.create_datatable()


    # you may want to include an existing datable to get the JS only once but 
    datatable = DataTable(df, 'data', datatable=datatable.datatable)
    datatable.datatable.set_tooltips_to_column('c','a')
    js = datatable.create_javascript_function()
    html = datatable.create_datatable()

