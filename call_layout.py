from app import app
from dash.dependencies import Input, Output

import pandas as pd
from collections import OrderedDict
import pandas as pd
import sqlite3
import dash_table as dt

from flask import Flask
import sqlite3
import dash
from dash import dash_table
import dash_bootstrap_components as dbc
from dash import dcc
from dash import html
import numpy as np
import plotly.graph_objs as go
from dash.dependencies import Input, Output


@app.callback(
    Output("tab_content", "children"), 
    [Input("tabs", "active_tab")])
def switch_tab(at):
    if at == "available_dbs":
        return tab_available_db_content 
    elif at == "annotations":
        return tab_annotation_content
    elif at == "sequences":
        return html.P("Seqs tab...")
    return html.P("This shouldn't ever be displayed...")

tab_annotation_content = html.Div([
    dcc.Textarea(
        id='gene-list',
        value='Textarea content initialized\nwith multiple lines of text',
        style={'width': '100%', 'height': 100},
        ),

    html.Div(id="annotation_table"),
    ])


def show_available_species():
    con = sqlite3.connect('SQNce.db')
    cursorObj = con.cursor()
    
    output_df = pd.read_sql_query("SELECT * FROM species", con)
    return dt.DataTable(
            columns=[{"name": i, "id": i} for i in output_df.columns],
            data=output_df.to_dict('records'))
tab_available_db_content = show_available_species()

#html.Div([
#    html.Div(id="available_dbs_table"),
#    ])


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
    Input('gene-list', 'value'))
def get_gene_list(value):
    #con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db')
    con = sqlite3.connect('SQNce.db')
    try:
        gene_list = value.split("\n")
        if len(gene_list) > 50:
            gene_list = gene_list[:50]
        if gene_list[-1]=="":
            gene_list = gene_list[:-1]
    except:
        return(html.P("Something did not work reading the gene list"))
    try:
        output_df = pd.DataFrame.from_dict(annotation_select(con, gene_list), orient="index").reset_index()
        output_df.columns = ["GeneID", "annotation"]
        return dt.DataTable(
            columns=[{"name": i, "id": i} for i in output_df.columns],
            data=output_df.to_dict('records'))
    except:
        return(html.P("Something did not work returning the table"))