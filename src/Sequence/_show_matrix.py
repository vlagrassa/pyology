def pretty_table_from_array(data_array, row_labels, col_labels):
    """Show an HTML table from a 2d numpy array"""
    from IPython.core.display import HTML,display
    import pandas as pd
    df = pd.DataFrame(data_array,index=row_labels,columns=col_labels)
    table_html = df.to_html()
    return HTML(table_html)

def display_table(data_array, row_labels, col_labels):
    display(pretty_table_from_array(data_array, row_labels, col_labels))