import sqlite3
import os 
import dash_bootstrap_components as dbc
from dash import Dash, dcc, html, Input, Output, State, dash_table


def register_callbacks(dashapp):
    import pandas as pd
    from collections import OrderedDict
    import dash_bootstrap_components as dbc
    from dash import Dash, dcc, html, Input, Output, State, dash_table, callback_context

    import csv

    from flask import Flask
    import sqlite3
    import numpy as np
    import plotly.graph_objs as go

    import os
    import zlib

    # Set the cwd and SQNce.db path to be platform independent  
    if "plantapp" not in os.getcwd().lower():
        cwd = "/home/eporetsky/plantapp" # for server hosting
    else:
        cwd = os.getcwd() # for personal computer
    sqnce_path = os.path.join(cwd, "SQNce.db")
    
    @dashapp.callback(Output('page-content', 'children'),
                      Input('sqnce_url', 'pathname'))
    def display_page(pathname):
        pathname = pathname.split("/")[-1]
        print("Current SQNce URL is:", pathname)
        if pathname == "available_dbs":
            return tab_avaialble_dbs
        elif pathname == "annotations":
            return tab_annotation_content
        elif pathname == "families_familyIDs":
            return tab_family_familyIDs_content
        elif pathname == "families_geneIDs":
            return tab_family_geneIDs_content
        elif pathname == "proteins":
            return tab_protein_content
        elif pathname == "promoters":
            return tab_promoter_content
        elif pathname == "coordinates":
            return tab_coordinates_content
        elif pathname == "omics":
            return tab_omics_content

    # Button trigger context: https://stackoverflow.com/questions/62119605/dash-how-to-callback-depending-on-which-button-is-being-clicked
    @dashapp.callback(
        Output("tab_content", "children"),
        [Input("available_dbs", "n_clicks"),
        Input("annotations", "n_clicks"),
        Input("families_familyIDs", "n_clicks"),
        Input("families_geneIDs", "n_clicks"),
        Input("coordinates", "n_clicks"),
        Input("proteins", "n_clicks"),
        Input("promoters", "n_clicks"),
        Input("omics", "n_clicks"),]
        )
    def switch_tab(*args):
        print()
        trigger = callback_context.triggered[0]
        trigger = trigger["prop_id"].split(".")[0]
        if trigger == "available_dbs":
            return tab_avaialble_dbs
        elif trigger == "annotations":
            return tab_annotation_content
        elif trigger == "families_familyIDs":
            return tab_family_familyIDs_content
        elif trigger == "families_geneIDs":
            return tab_family_geneIDs_content
        elif trigger == "proteins":
            return tab_protein_content
        elif trigger == "promoters":
            return tab_promoter_content
        elif trigger == "coordinates":
            return tab_coordinates_content
        elif trigger == "omics":
            return tab_omics_content
        return html.P("This shouldn't ever be displayed...")

    ###############################################################################
    #                                Helper Functions
    ###############################################################################
    def distinct_db_vals(db, table, column, custom_vals=[], return_ls=False):
        # Input is the column to select and from which table
        # Returns a list of all values in a specific table from SQNce.db
        # Custom vals are added to the front using nested list of [label, value]
        # Use return_ls for clean list of distinct values, not for dropdown menu
        ls = [{ 'label': label, 'value': val} for label, val in custom_vals]
        con = sqlite3.connect(db) # deploy with this
        cursorObj = con.cursor()
        distinct_df = pd.read_sql_query('''SELECT DISTINCT {0} 
                                        FROM {1}'''.format(column, table), con)
        if return_ls:
            return(distinct_df[column].to_list())
        for name in distinct_df[column]:
            ls.append({'label': name, 'value': name})
        return(ls)

    ###############################################################################
    #                                Coordinates
    ###############################################################################

    tab_coordinates_content = html.Div([

        # Get a list of unique genotypes from the db
        dcc.Dropdown(
            id='coordinates_genotypes_dropdown',
            options=distinct_db_vals(sqnce_path, "gene_coordinates", "genotype_id", 
                                     [["None", "None"]]),
            value="None"
        ),
        
        # Select how many bp from each side to search
        dcc.Dropdown(
            id='coordinates_range_dropdown',
            options=[
                {'label': '1,000bp', 'value': 1000},
                {'label': '10,000bp', 'value': 10000},
                {'label': '100,000bp', 'value': 100000}
            ],
            value=1000
        ),

        dcc.Dropdown(
            id='coordinates_dedupe_dropdown',
            options=[
                {'label': "Remove duplicates", 'value': True},
                {'label': "Keep all values", 'value': False}
            ],
            value="None"
        ),

        dcc.Textarea(
            id='coordinates_list',
            value='Textarea content initialized\nwith multiple lines of text',
            style={'width': '100%', 'height': 100},
            ),

        dbc.Row(dbc.Col(
            [
            dbc.Button("Download coordinate table", color="primary", id="btn_download_coordinate_table", className="mr-1"),
            dcc.Download(id="download_coordinate_table"),
        ])),
        dcc.Clipboard(id="coordinate_table_copy", style={"fontSize":20}),

        html.Div(id="coordinates_table"),
    ]),

    # Query to find neighboring genes
    def get_SNP_neighbors(genotype, chromsome, coordinate, distance):
        con = sqlite3.connect(sqnce_path) # deploy with this
        cursorObj = con.cursor()
        df = pd.read_sql_query('''SELECT * 
                        FROM gene_coordinates 
                        WHERE genotype_id = "{0}"
                        AND gene_chr = "{1}"
                        AND gene_start BETWEEN {2} AND {3}
                        
                        UNION ALL
                        
                        SELECT * 
                        FROM gene_coordinates 
                        WHERE genotype_id = "{0}"
                        AND gene_chr = "{1}"
                        AND gene_end BETWEEN {2} AND {3}
                        '''.format(genotype, chromsome, coordinate-distance, coordinate+distance), con)
        # Should check why it returns the same row twice, probably need to correct the query
        df = df.drop_duplicates()
        df.insert(0, 'Query', pd.Series(["_".join([chromsome, str(coordinate)]) for x in range(len(df.index))]))
        return(df)

    @dashapp.callback(
        Output('coordinates_table', 'children'),
        Input('coordinates_list', 'value'),
        Input("coordinates_genotypes_dropdown", 'value'),        
        Input("coordinates_range_dropdown", 'value'),
        Input("coordinates_dedupe_dropdown", 'value'),
    )
    def get_coordinates_gene_list(value, genotype, distance, dedupe):
        # Dedupe - remove duplicates that might appear of entries are close-by
        con = sqlite3.connect(sqnce_path) # deploy with this
        cursorObj = con.cursor()
        try:
            rows = value.split("\n")
            coordinate_list = [row.split("\t") for row in rows]
            if coordinate_list[-1]==[""]:
                coordinate_list = coordinate_list[:-1]
            for row in coordinate_list:
                if len(row) != 2:
                    return(html.P("Number of columns is not 2. Use tab-seperated values."))
        except:
            return(html.P("Something did not work while reading the coordinate list"))
        try:
            df = pd.DataFrame()
            # Check the neighboring genes for each row in the input data
            for row in coordinate_list:
                df = pd.concat([df, get_SNP_neighbors(genotype, str(row[0]), int(row[1]), distance)])
            # Since some genes can overlap, remove duplicated values and eep first
            if dedupe:
                df = df.drop_duplicates(subset=["gene_id"]).reset_index()

            # Add a column of gene annotations when available
            df["annotation"] = annotation_select(con, df["gene_id"].to_list())

            # Export the dast datatable to the content container
            return dash_table.DataTable(
                id="coordinate_table_state",
                columns=[{"name": i, "id": i} for i in df.columns],
                data=df.to_dict('records'),
                style_cell={'textAlign': 'left',
                            'overflow': 'hidden',
                            'textOverflow': 'ellipsis',
                            'maxWidth': 0,},
                editable=True,
                row_deletable=True,)
        except:
            return(html.P("Something did not work returning the table"))


    @dashapp.callback(
        Output("coordinate_table_copy", "content"),
        Input("coordinate_table_copy", "n_clicks"),
        State("coordinate_table_state", "data"),
    )
    def coordinates_copy(_, data):
        df = pd.DataFrame(data)
        # See options for .to_csv() or .to_excel() or .to_string() in the  pandas documentation
        return df.to_csv(index=False, sep="\t", quoting=csv.QUOTE_NONE)  # includes headers

    @dashapp.callback(
        Output("download_coordinate_table", "data"),
        Input("btn_download_coordinate_table", "n_clicks"),
        State("coordinate_table_state", "data"),
        prevent_initial_call=True,
    )
    def download_coordinate_fasta(n_clicks, data):
        dff = pd.DataFrame(data)
        output = dff.to_csv(index=False, header=True, sep="\t", quoting=csv.QUOTE_NONE)
        return dict(content=output, filename="coordinate_output.tsv")

    ###############################################################################
    #                            Analyzed Studies
    ###############################################################################

    def get_studies_in_species(selected_species):
        con = sqlite3.connect(sqnce_path) # deploy with this
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

    @dashapp.callback(
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
        con = sqlite3.connect(sqnce_path) # deploy with this
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

    @dashapp.callback(
        Output('samples_in_selected_study', 'children'),
        [Input('studies-species-dropdown', 'value'),
        Input('studies-in-species-dropdown', 'value')])
    def samples_in_selected_study(species, study):
        if study=="None":
            return(html.P("No study is selected."))
        else: 
            samples_df = get_fastq_table(species, study)
            return dash_table.DataTable(
                columns=[{"name": i, "id": i} for i in samples_df.columns],
                data=samples_df.to_dict('records'))

    tab_omics_content = html.Div([
        dcc.Dropdown(
            id='studies-species-dropdown',
            options=distinct_db_vals(sqnce_path, "studies", "scientific_name", 
                                     [["None", "None"]]),
            value="None"
            ),
        html.Div(id="studies-in-species-dropdown-div"),
        html.Div(id="samples_in_selected_study"),
        ])

    ###############################################################################
    #                            Availalbe Datasets
    ###############################################################################

    tab_avaialble_dbs = html.Div([
        # Get a list of unique genotypes from the db
        dcc.Dropdown(
            id='avaialble_dbs_dropdown',
            options=[{'label': fl.split(".")[0], 'value': fl.split(".")[0]} for fl in os.listdir(os.path.join(cwd,"init"))],
            value="None"
        ),
        html.Div(id="available_dbs_table"),
    ])

    @dashapp.callback(
        Output("available_dbs_table", "children"),
        Input("avaialble_dbs_dropdown", "value"),
    )
    def select_available_dbs(value):
        try:
            df = pd.read_csv(os.path.join(cwd,"init/",value+".tsv"), sep="\t")
            return dash_table.DataTable(
                id="available_dbs_table_state",
                columns=[{"name": i, "id": i} for i in df.columns],
                data=df.to_dict('records'),
                style_cell={'textAlign': 'left',
                            'overflow': 'hidden',
                            'textOverflow': 'ellipsis',
                            'maxWidth': 0,},
                editable=True,
                row_deletable=True,)
        except:
            return(html.P("Something did not work returning the table"))
    
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
        ls = []
        for entity in entity_list:
            cursorObj = con.cursor()
            cursorObj.execute('''SELECT gene_id, gene_annotation 
                                FROM gene_annotations 
                                WHERE gene_id =  ?  ''', (entity,))
            # (name,) - need the comma to treat it as a single item and not list of letters
            selected = cursorObj.fetchall()
            if selected == []:
                ls.append("Gene not found")
            else:
                ls.append(selected[0][1])    
        return(ls)

    @dashapp.callback(
        Output('annotation_table', 'children'),
        Input('gene-list', 'value'))
    def get_gene_list(value):
        con = sqlite3.connect(sqnce_path) # deploy with this
        try:
            gene_list = value.split("\n")
            if len(gene_list) > 500:
                gene_list = gene_list[:500]
            if gene_list[-1]=="":
                gene_list = gene_list[:-1]
        except:
            return(html.P("Something did not work reading the gene list"))
        try:
            output_df = pd.DataFrame({"GeneID": gene_list, "annotation": annotation_select(con, gene_list) })
            output_df.columns = ["GeneID", "annotation"]
            return dash_table.DataTable(
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
    #                                Family Annotations
    ###############################################################################

    tab_family_familyIDs_content = html.Div([
        dcc.Dropdown(
            id='family_select_species_dropdown',
            options=distinct_db_vals(sqnce_path, "gene_families", "species_id", 
                                    [["All", "all"]]),
            value=None,
            multi=True,
            searchable=True
        ),
        dcc.Dropdown(
            id='family_select_name_dropdown',
            options=distinct_db_vals(sqnce_path, "gene_families", "family_name"),
            value=None,
            multi=True,
            searchable=True
        ),
        dbc.Button("Get selected gene families", id='family_get_selected_btn', 
                    outline=True, color="primary", className="me-1"),
                
        html.Div(id="family_familyIDs_table"),
        ])
    
    @dashapp.callback(
        Output("family_familyIDs_table", "children"),
        Input("family_get_selected_btn", "n_clicks"),
        State("family_select_species_dropdown", "value"),
        State("family_select_name_dropdown", "value"),
    )
    def family_familyIDs_table(_, genotypes, families):
        con = sqlite3.connect(sqnce_path) # deploy with this
        cursorObj = con.cursor()
        if  "all" in genotypes:
            genotypes = distinct_db_vals(sqnce_path, "gene_families", "species_id", return_ls=True) 
        genotypes = str("','".join(genotypes))
        families = str("','".join(families))
        df = pd.read_sql_query("""SELECT protein_id, family_name 
                                FROM gene_families
                                WHERE species_id IN ('{0}') AND family_name IN ('{1}')""".format(genotypes, families), con)
        try:
            return dash_table.DataTable(
                id="family_table_state",
                columns=[{"name": i, "id": i} for i in df.columns],
                data=df.to_dict('records'),
                style_cell={'textAlign': 'left',
                            'overflow': 'hidden',
                            'maxWidth': 0,},
                editable=True,
                row_deletable=True,)
        except:
            return(html.P("Something did not work returning the genotype-family table"))
    
    #####################################################################################

    tab_family_geneIDs_content = html.Div([
        dcc.Textarea(
            id='family_gene_list',
            value='Insert gene or gene-family list here, one id per row.',
            style={'width': '100%', 'height': 100},
            ),
            html.Div(id="family_geneIDs_table"),
        ])
    
    def family_gene_select(con, gene_list):
        # Use an input list of genes to find their family assignments
        ls = []
        for gene in gene_list:
            cursorObj = con.cursor()
            cursorObj.execute('''SELECT protein_id, family_name 
                                FROM gene_families
                                WHERE protein_id =  ? ''', (gene,))
            selected = cursorObj.fetchall()
            if selected == []:
                ls.append("Gene not found")
            else:
                ls.append(selected[0][1])    
        return(ls)
    
    @dashapp.callback(
        Output("family_geneIDs_table", "children"),
        Input('family_gene_list', 'value'),
    )
    def family_geneIDs_table(value):
        con = sqlite3.connect(sqnce_path) # deploy with this
        try:
            gene_list = value.split("\n")
            if len(gene_list) > 1000:
                gene_list = gene_list[:1000]
            if gene_list[-1]=="":
                gene_list = gene_list[:-1]
        except:
            return(html.P("Something did not work reading the gene list"))
        try:
            output_df = pd.DataFrame({"GeneID": gene_list, "FamilyID": family_gene_select(con, gene_list) })
            # output_df.columns = ["GeneID", "FamilyID"]
            return dash_table.DataTable(
                id="family_table_state",
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
            #print(od)
        return(od)

    @dashapp.callback(
        Output('protein_seq_table', 'children'),
        Input('protein-gene-list', 'value'))
    def get_protein_list(value):
        con = sqlite3.connect(sqnce_path) # deploy with this
        try:
            gene_list = value.split("\n")
            if len(gene_list) > 500:
                gene_list = gene_list[:500]
            if gene_list[-1]=="":
                gene_list = gene_list[:-1]
        except:
            return(html.P("Something did not work reading the gene list"))
        try:
            output_df = pd.DataFrame.from_dict(proteins_select(con, gene_list), orient="index").reset_index()
            output_df.columns = ["GeneID", "Sequence"]
            return dash_table.DataTable(
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
        
    @dashapp.callback(
        Output("protein_table_copy", "content"),
        Input("protein_table_copy", "n_clicks"),
        State("protein_table_state", "data"),
    )
    def custom_copy(_, data):
        dff = pd.DataFrame(data)
        # See options for .to_csv() or .to_excel() or .to_string() in the  pandas documentation
        return dff.to_csv(index=False)  # includes headers

    @dashapp.callback(
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

    @dashapp.callback(
        Output('promoter_seq_table', 'children'),
        Input('promoter-gene-list', 'value'),
        Input('promoters-kind-dropdown', 'value'))
    def get_promoter_list(gene_list, promoter_kind):
        con = sqlite3.connect(sqnce_path) # deploy with this
        try:
            gene_list = gene_list.split("\n")
            if len(gene_list) > 10000:
                gene_list = gene_list[:10000]
            if gene_list[-1]=="":
                gene_list = gene_list[:-1]
        except:
            return(html.P("Something did not work reading the gene list +"+str(gene_list)+promoter_kind))
        try:
            output_df = pd.DataFrame.from_dict(promoter_select(con, gene_list, promoter_kind), orient="index").reset_index()
            output_df.columns = ["GeneID", "Sequence"]
            return dash_table.DataTable(
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
        
    @dashapp.callback(
        Output("promoter_table_copy", "content"),
        Input("promoter_table_copy", "n_clicks"),
        State("promoter_table_state", "data"),
    )
    def custom_promoter_copy(_, data):
        df = pd.DataFrame(data)
        # See options for .to_csv() or .to_excel() or .to_string() in the  pandas documentation
        return df.to_csv(index=False)  # includes headers

    @dashapp.callback(
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

    
        