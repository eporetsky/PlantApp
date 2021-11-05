from app import app
from dash.dependencies import Input, Output

import pandas as pd
from collections import OrderedDict
import pandas as pd
import sqlite3
from dash import dash_table as dt
from dash import dash_table

import csv

from flask import Flask
import sqlite3
import dash
import dash_bootstrap_components as dbc
from dash import Dash, dcc, html, Input, Output, State, dash_table

import numpy as np
import plotly.graph_objs as go
from dash.dependencies import Input, Output

import os
import zlib

@app.callback(
    Output("tab_content", "children"), 
    [Input("tabs", "active_tab")])
def switch_tab(at):
    if at == "available_dbs":
        return show_available_species() 
    elif at == "annotations":
        return tab_annotation_content
    elif at == "proteins":
        return tab_protein_content
    elif at == "promoters":
        return tab_promoter_content
    elif at == "omics":
        return tab_omics_content
    return html.P("This shouldn't ever be displayed...")

###############################################################################
#                            Analyzed Studies
###############################################################################

def get_studies_species():
    # Returns a list of species names of all analyzed studies from SQNce.db
    #con = sqlite3.connect(os.getcwd()+'/SQNce.db')
    con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db')
    cursorObj = con.cursor()
    distinct_species_df = pd.read_sql("""
                SELECT DISTINCT studies.scientific_name
                FROM studies
                """, con = con)
    dropdown_species_list = [{'label': "None", 'value': "None"}]
    for name in distinct_species_df["scientific_name"]:
        dropdown_species_list.append({'label': name, 'value': name})
    return(dropdown_species_list)

def get_studies_in_species(selected_species):
    #con = sqlite3.connect(os.getcwd()+'/SQNce.db')
    con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db')
    cursorObj = con.cursor()
    sql_query = """SELECT studies.study_accession
                FROM studies
                WHERE studies.scientific_name='%s'
    """%(selected_species)
    studies_in_species_df = pd.read_sql(sql_query, con = con)
    studies_in_species_list = [{'label': "Select all", 'value': "all"}]
    for name in studies_in_species_df["study_accession"]:
        studies_in_species_list.append({'label': name, 'value': name})
    return(studies_in_species_list)

@app.callback(
    Output('studies-in-species-dropdown-div', 'children'),
    Input('studies-species-dropdown', 'value'))
def studies_in_species(value):
    if value == "None":
        return(html.P("No species is selected."))
    else:
        return(dcc.Dropdown(
            id='studies-in-species-dropdown',
            options=get_studies_in_species(value),
            value="None"
            ),
        )
        
def get_fastq_table(selected_species, selected_studies):
    #con = sqlite3.connect(os.getcwd()+'/SQNce.db')
    con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db')
    cursorObj = con.cursor()
    if selected_studies=="all":
        query_line = "studies.scientific_name='{0}'".format(selected_species)
    else:
        query_line = "fastq.study_accession='{0}'".format(selected_studies)
    sql_query = """SELECT *
                FROM studies, fastq
                WHERE studies.study_accession=fastq.study_accession and %s
                """%(query_line)
    return(pd.read_sql_query(sql_query, con))

@app.callback(
    Output('samples_in_selected_study', 'children'),
    [Input('studies-species-dropdown', 'value'),
    Input('studies-in-species-dropdown', 'value')])
def samples_in_selected_study(species, study):
    if study=="None":
        return(html.P("No study is selected."))
    else: 
        samples_df = get_fastq_table(species, study)
        return dt.DataTable(
            columns=[{"name": i, "id": i} for i in samples_df.columns],
            data=samples_df.to_dict('records'))

tab_omics_content = html.Div([
    dcc.Dropdown(
        id='studies-species-dropdown',
        options=get_studies_species(),
        value="None"
        ),
    html.Div(id="studies-in-species-dropdown-div"),
    html.Div(id="samples_in_selected_study"),
    ])


###############################################################################
#                            Availalbe Datasets
###############################################################################

def show_available_species():
    #con = sqlite3.connect(os.getcwd()+'/SQNce.db')
    con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db')
    cursorObj = con.cursor()
    
    output_df = pd.read_sql_query("SELECT * FROM species", con)
    return dt.DataTable(
            columns=[{"name": i, "id": i} for i in output_df.columns],
            data=output_df.to_dict('records'))

###############################################################################
#                                Annotations
###############################################################################

tab_annotation_content = html.Div([
    dcc.Textarea(
        id='gene-list',
        value='Textarea content initialized\nwith multiple lines of text',
        style={'width': '100%', 'height': 100},
        ),

    html.Div(id="annotation_table"),
    ])

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
    con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db')
    #con = sqlite3.connect(os.getcwd()+'/SQNce.db')
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
            data=output_df.to_dict('records'),
            style_cell={'textAlign': 'left',
                        'overflow': 'hidden',
                        'textOverflow': 'ellipsis',
                        'maxWidth': 0,},
            editable=True,
            row_deletable=True,)
    except:
        return(html.P("Something did not work returning the table"))
    
###############################################################################
#                                Protein Sequences
###############################################################################

tab_protein_content = html.Div([
    dcc.Textarea(
        id='protein-gene-list',
        value='Textarea content initialized\nwith multiple lines of text',
        style={'width': '100%', 'height': 100},
        ),
    dbc.Row(dbc.Col(
        [
         dbc.Button("Download fasta with gene IDs", color="primary", id="btn_download_fasta_geneIDs", className="mr-1"),
         dcc.Download(id="download_fasta_geneIDs"),
         dbc.Button("Download fasta with symbols", color="primary", className="mr-1"),
    ])),
    dcc.Clipboard(id="protein_table_copy", style={"fontSize":20}),
    html.Div(id="protein_seq_table"),
    ])

def proteins_select(con, entity_list):
    od = OrderedDict()
    for entity in entity_list:
        cursorObj = con.cursor()
        cursorObj.execute('''SELECT protein_id, protein_sequence
                             FROM protein_seqs
                             WHERE protein_id =  ?  ''', (entity,))
        # (name,) - need the comma to treat it as a single item and not list of letters
        selected = cursorObj.fetchall()[0]
        od[selected[0]] = zlib.decompress(selected[1]).decode('utf-8')[:-1] 
        #od[selected[0]] = selected[1][:-1]
        print(od)
    return(od)

@app.callback(
    Output('protein_seq_table', 'children'),
    Input('protein-gene-list', 'value'))
def get_protein_list(value):
    con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db')
    #con = sqlite3.connect(os.getcwd()+'/SQNce.db')
    try:
        gene_list = value.split("\n")
        if len(gene_list) > 50:
            gene_list = gene_list[:50]
        if gene_list[-1]=="":
            gene_list = gene_list[:-1]
    except:
        return(html.P("Something did not work reading the gene list"))
    try:
        output_df = pd.DataFrame.from_dict(proteins_select(con, gene_list), orient="index").reset_index()
        output_df.columns = ["GeneID", "Sequence"]
        return dt.DataTable(
            id="protein_table_state",
            columns=[{"name": i, "id": i} for i in output_df.columns],
            data=output_df.to_dict('records'),
            style_cell={'textAlign': 'left',
                        'overflow': 'hidden',
                        'maxWidth': 0,},
            editable=True,
            row_deletable=True,)
    except:
        return(html.P("Something did not work returning the table"))
    
@app.callback(
    Output("protein_table_copy", "content"),
    Input("protein_table_copy", "n_clicks"),
    State("protein_table_state", "data"),
)
def custom_copy(_, data):
    dff = pd.DataFrame(data)
    # See options for .to_csv() or .to_excel() or .to_string() in the  pandas documentation
    return dff.to_csv(index=False)  # includes headers

@app.callback(
    Output("download_fasta_geneIDs", "data"),
    Input("btn_download_fasta_geneIDs", "n_clicks"),
    State("protein_table_state", "data"),
    prevent_initial_call=True,
)
def func(n_clicks, data):
    dff = pd.DataFrame(data)
    dff ["GeneID"] = ">"+dff ["GeneID"]
    dff = dff.agg('\n'.join, axis=1) 
    # By adding \n when aggregating the two columns we bascially create a fasta file
    # which doesn't allow \n in it so have to use escapechar: https://github.com/pandas-dev/pandas/issues/16298
    output = dff.to_csv(index=False, header=False, escapechar="#", quoting=csv.QUOTE_NONE, line_terminator="\n").replace("#", "")
    return dict(content=output, filename="output.fasta")


###############################################################################
#                                Promoter Sequences
###############################################################################

tab_promoter_content = html.Div([
    dcc.Textarea(
        id='promoter-gene-list',
        value='Insert list of genes to\nget a list of promoter sequences',
        style={'width': '100%', 'height': 100},
        ),
    dbc.Row(dbc.Col(
        [
         dbc.Button("Download fasta with gene IDs", color="primary", id="btn_promoter_download_fasta_geneIDs", className="mr-1"),
         dcc.Download(id="download_promoter_fasta_geneIDs"),
         dbc.Button("Download fasta with symbols", color="primary", className="mr-1"),
    ])),
    dcc.Clipboard(id="promoter_table_copy", style={"fontSize":20}),
    html.Div(id="promoter_seq_table"),
    ])

def promoter_select(con, entity_list):
    od = OrderedDict()
    for entity in entity_list:
        cursorObj = con.cursor()
        cursorObj.execute('''SELECT protein_id, promoter_sequence
                             FROM promoter_seqs
                             WHERE protein_id =  ?  ''', (entity,))
        # (name,) - need the comma to treat it as a single item and not list of letters
        selected = cursorObj.fetchall()[0]
        od[selected[0]] = zlib.decompress(selected[1]).decode('utf-8')[:-1]
        #od[selected[0]] = selected[1][:-1]
    return(od)

@app.callback(
    Output('promoter_seq_table', 'children'),
    Input('promoter-gene-list', 'value'))
def get_promoter_list(value):
    con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db')
    #con = sqlite3.connect(os.getcwd()+'/SQNce.db')
    try:
        gene_list = value.split("\n")
        if len(gene_list) > 1000:
            gene_list = gene_list[:1000]
        if gene_list[-1]=="":
            gene_list = gene_list[:-1]
    except:
        return(html.P("Something did not work reading the gene list"))
    try:
        output_df = pd.DataFrame.from_dict(promoter_select(con, gene_list), orient="index").reset_index()
        output_df.columns = ["GeneID", "Sequence"]
        return dt.DataTable(
            id="promoter_table_state",
            columns=[{"name": i, "id": i} for i in output_df.columns],
            data=output_df.to_dict('records'),
            style_cell={'textAlign': 'left',
                        'overflow': 'hidden',
                        'maxWidth': 0,},
            editable=True,
            row_deletable=True,)
    except:
        return(html.P("Something did not work returning the table"))
    
@app.callback(
    Output("promoter_table_copy", "content"),
    Input("promoter_table_copy", "n_clicks"),
    State("promoter_table_state", "data"),
)
def custom_promoter_copy(_, data):
    df = pd.DataFrame(data)
    # See options for .to_csv() or .to_excel() or .to_string() in the  pandas documentation
    return df.to_csv(index=False)  # includes headers

@app.callback(
    Output("download_promoter_fasta_geneIDs", "data"),
    Input("btn_promoter_download_fasta_geneIDs", "n_clicks"),
    State("promoter_table_state", "data"),
    prevent_initial_call=True,
)
def download_promoter_fasta(n_clicks, data):
    dff = pd.DataFrame(data)
    dff ["GeneID"] = ">"+dff ["GeneID"]
    dff = dff.agg('\n'.join, axis=1) 
    # By adding \n when aggregating the two columns we bascially create a fasta file
    # which doesn't allow \n in it so have to use escapechar: https://github.com/pandas-dev/pandas/issues/16298
    output = dff.to_csv(index=False, header=False, escapechar="#", quoting=csv.QUOTE_NONE, line_terminator="\n").replace("#", "")
    return dict(content=output, filename="output.fasta")