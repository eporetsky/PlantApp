import sqlite3
import os 
con = sqlite3.connect(os.getcwd()+'/SQNce.db')

def register_callbacks(dashapp):
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

    #con = sqlite3.connect(os.getcwd()+'/SQNce.db')
    #con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db')

    ###############################################################################
    #                            Analyzed Studies
    ###############################################################################

    def get_studies_species():
        # Returns a list of species names of all analyzed studies from SQNce.db
        #con = sqlite3.connect(os.getcwd()+'/SQNce.db')
        #con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db')
        global con
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
        #con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db')
        global con
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
            
    def get_fastq_table(selected_species, selected_studies):
        #con = sqlite3.connect(os.getcwd()+'/SQNce.db')
        #con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db')
        global con
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
        #con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db')
        global con
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

    ###############################################################################
    #                                Promoter Sequences
    ###############################################################################

    tab_promoter_content = html.Div([
        dcc.Dropdown(
            id='promoters-kind-dropdown',
            options=[{'label': "TSS", 'value': "TSS"}, {'label': "ATG", 'value': "ATG"}],
            value="TSS"
            ),
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

    def promoter_select(con, entity_list, promoter_kind):
        od = OrderedDict()
        for entity in entity_list:
            cursorObj = con.cursor()
            cursorObj.execute('''SELECT protein_id, promoter_sequence
                                FROM promoter_seqs
                                WHERE protein_id =  ? and promoter_kind = ?''', [entity, promoter_kind])
            # (name,) - need the comma to treat it as a single item and not list of letters
            selected = cursorObj.fetchall()[0]
            od[selected[0]] = zlib.decompress(selected[1]).decode('utf-8')[:-1]
            #od[selected[0]] = selected[1][:-1]
        return(od)