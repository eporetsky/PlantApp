from app import app
import sqlite3
import dash_table as dt
import pandas as pd
from collections import OrderedDict
from dash.dependencies import Input, Output


def annotation_select(con, entity_list):
    od = OrderedDict()
    for entity in entity_list:
        cursorObj = con.cursor()
        cursorObj.execute('''SELECT gene_id, gene_annotation 
                             FROM gene_annotations 
                             WHERE gene_id =  ?  ''', (entity,))
        # (name,) - need the comma to treat it as a single item and not list of letters
        selected = cursorObj.fetchall()[0]
        od[selected[0]] = selected[1]
    return(od)



@app.callback(
    Output('annotation_table', 'children'),
    Input('gene-list', 'value')
)
def get_gene_list(value):
	con = sqlite3.connect('SQNce.db')
	try:
		gene_list = value.split("\n")
	except:
		gene_list = ["Zm00001d027231"]
	output_df = pd.DataFrame.from_dict(annotation_select(con, gene_list), orient="index").reset_index()
	output_df.columns = ["GeneID", "annotation"]
	
	return dt.DataTable(
    		columns=[{"name": i, "id": i} for i in output_df.columns],
    		data=output_df.to_dict('records'),
    		)