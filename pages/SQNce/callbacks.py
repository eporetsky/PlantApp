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
    import duckdb
    #import mysql.connector

    # Set the cwd and SQNce.db path to be platform independent  
    if "plantapp" not in os.getcwd().lower():
        cwd = "/home/eporetsky/plantapp" # for server hosting
    else:
        cwd = os.getcwd() # for personal computer
    sqnce_path = os.path.join(cwd, "SQNce.db")
    
    #con = mysql.connector.connect(
    #                host='username.mysql.pythonanywhere-services.com', 
    #                db='username$DatabaseName', 
    #                user='username', 
    #                password='password'
    #                ) 


    @dashapp.callback(Output('tab_content', 'children'),
                      Input('url', 'pathname'))
    def display_page_content(pathname):
        #print(@app.context_processor)
        pathname = pathname.split("/")[-1]
        if pathname == "available_dbs":
            return tab_avaialble_dbs
        elif pathname == "annotations":
            return tab_annotation_content
        elif pathname == "symbols":
            return tab_symbols_content
        elif pathname == "families_familyIDs":
            return tab_family_familyIDs_content
        elif pathname == "families_geneIDs":
            return tab_family_geneIDs_content
        elif pathname == "orthogroups":
            return tab_orthogroups_content        
        elif pathname == "BBHs":
            return tab_BBHs_content
        elif pathname == "proteins":
            return tab_protein_content
        elif pathname == "promoters":
            return tab_promoter_content
        elif pathname == "gene_coordinates":
            return tab_gene_coordiantes_content
        elif pathname == "genes_near_coordinates":
            return tab_genes_near_coordinates_content
        elif pathname == "transcriptomics":
            return tab_transcriptomics_content
        elif pathname == "omics":
            return tab_omics_content

    ###############################################################################
    #                                Helper Functions
    ###############################################################################
    def distinct_db_vals(db, table, column, custom_vals=[], return_ls=False):
        # Input is the column to select and from which table
        # Returns a list of all values in a specific table from SQNce.db
        # Custom vals are added to the front using nested list of [label, value]
        # Use return_ls for clean list of distinct values, not for dropdown menu
        try:
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
        except:
            return None

    ###############################################################################
    #                               Gene Coordinates
    ###############################################################################

    tab_gene_coordiantes_content = html.Div([
        html.Div(
            className="row", children=[
                html.Div(children=[html.P("Add gene IDs below to get their genomic coordiantes:", style={'marginTop': 10, 'marginBottom': 10})], 
                    className='six columns', style=dict(width='25%')), 
            ], style=dict(display='flex')),
        
        html.Div(
            className="row", children=[
                html.Div(children=[dcc.Textarea(id='gene_coordinates_gene_list',
                    value='Add multiple gene IDs below, 1 gene ID per row\n(maximum gene IDs per query is 60,000)',
                    style={'width': '80%', 'height': 60, 'marginBottom': 10, 'marginTop': 0},
                    )], className='six columns', style=dict(width='40%')),
            ], style=dict(display='flex')),

        html.Div(
            className="row", children=[
                html.Div(children=[dbc.Button("Find gene coordiantes", 
                                    color="primary", id="btn_submit_gene_coordinates_gene_list", className="mr-1",
                                    style={'marginTop': 0, 'marginBottom': 20})], 
                    className='six columns', style=dict(width='25%')), 
            ], style=dict(display='flex')),

        html.Div(
            className="row", children=[
                    html.Div(children=[html.P("Copy table to clipboard:")],
                        className='six columns', style=dict(width='12%')),
                    html.Div(children=[dcc.Clipboard(id="gene_coordinates_copy", style={"fontSize":20})], 
                        className='six columns', style=dict(width='5%')),
            ], style=dict(display='flex')),

                
        html.Div(id="gene_coordinates_table"),
        ], style=dict(marginLeft=20, marginRight=20))

    @dashapp.callback(
        Output("gene_coordinates_copy", "content"),
        Input("gene_coordinates_copy", "n_clicks"),
        State("state_gene_coordinates_table", "data"),
        prevent_initial_call=True,
    )
    def orthogroup_genes_copy(_, df):
        df = pd.DataFrame(df)
        # See options for .to_csv() or .to_excel() or .to_string() in the  pandas documentation
        return df.to_csv(index=False)  # includes headers

    @dashapp.callback(
        Output('gene_coordinates_table', 'children'),
        Input('btn_submit_gene_coordinates_gene_list', 'n_clicks'),
        State('gene_coordinates_gene_list', 'value'),
        prevent_initial_call=True,)
    def get_gene_coordinates(_, gene_list):
        print("Getting gene coordiantes")
        con = sqlite3.connect(sqnce_path) # deploy with this
        try:
            gene_list = gene_list.split("\n")
            if len(gene_list) > 70000:
                gene_list = gene_list[:70000]
            if gene_list[-1]=="":
                gene_list = gene_list[:-1]

            
            gene_coordinates_list = []
            for gene_id in gene_list:
                cursorObj = con.cursor()
                cursorObj.execute('''SELECT gene_id, genome_id, gene_chr, gene_start, gene_end, gene_orientation
                                    FROM gene_coordinates
                                    WHERE gene_id =  ?  ''', (gene_id,))
                # (name,) - need the comma to treat it as a single item and not list of letters
                selected = cursorObj.fetchall()
                if selected == []:
                    continue # Skip if cannot find gene
                    #od[selected[0]] =  ...
                selected = selected[0] # Otherwise, get the value of the first returned row
                gene_coordinates_list.append(selected)
        except:
            return(html.P("Something did not work reading the gene list"))

        try:
            df = pd.DataFrame(gene_coordinates_list, columns=["gene_id", "genome_id", "gene_chr", "gene_start", "gene_end", "gene_orientation"])
            gene_coordiantes_table = dash_table.DataTable(
                id="state_gene_coordinates_table",
                columns=[{"name": i, "id": i} for i in df.columns],
                data=df.to_dict('records'),
                style_cell={'textAlign': 'left',
                            'overflow': 'hidden',
                            'maxWidth': 0,},
                sort_action='native',
                sort_mode='multi',
                editable=True,
                row_deletable=True,)

            return(gene_coordiantes_table)
        except:
            return(html.P("Something did not work reading the gene list"))

    ###############################################################################
    #                      Genes Near Coordinates
    ###############################################################################

    tab_genes_near_coordinates_content = html.Div([

        # Get a list of unique genotypes from the db
        dcc.Dropdown(
            id='coordinates_genotypes_dropdown',
            options=distinct_db_vals(sqnce_path, "gene_coordinates", "genome_id", 
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
        dbc.Button("Show example coordinates", color="primary", id="btn_coordinate_example", className="mr-1"),

        html.Div(id="coordinates_table"),
    ]),

    @dashapp.callback(
        Output('coordinates_genotypes_dropdown', 'value'),
        Output('coordinates_range_dropdown', 'value'),
        Output('coordinates_dedupe_dropdown', 'value'),
        Output('coordinates_list', 'value'),
        Input('btn_coordinate_example', 'n_clicks'),
        prevent_initial_call=True,)
    def symbolcoordinate_example(value):
        return "B73v4", 1000, True, '1	8574461'


    # Query to find neighboring genes
    def get_SNP_neighbors(genotype, chromsome, coordinate, distance):
        con = sqlite3.connect(sqnce_path) # deploy with this
        cursorObj = con.cursor()
        df = pd.read_sql_query('''SELECT * 
                        FROM gene_coordinates 
                        WHERE genome_id = "{0}"
                        AND gene_chr = "{1}"
                        AND gene_start BETWEEN {2} AND {3}
                        
                        UNION ALL
                        
                        SELECT * 
                        FROM gene_coordinates 
                        WHERE genome_id = "{0}"
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
        try:
            sql_query = """SELECT studies.study_accession
                    FROM studies
                    WHERE studies.scientific_name='%s'
            """%(selected_species)
            studies_in_species_df = pd.read_sql(sql_query, con = con)
            studies_in_species_list = [{'label': "Select all", 'value': "all"}]
            for name in studies_in_species_df["study_accession"]:
                studies_in_species_list.append({'label': name, 'value': name})
            return(studies_in_species_list)
        except:
            return None
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
            id='annotation_gene_list',
            value='Textarea content initialized\nwith multiple lines of text',
            style={'width': '100%', 'height': 100},
            ),

        dbc.Button("Show example annotaion", color="primary", id="btn_annotation_example", className="mr-1"),

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
        Output('annotation_gene_list', 'value'),
        Input('btn_annotation_example', 'n_clicks'),
        prevent_initial_call=True,)
    def annotation_example(value):
        return("Zm00001d021929\nZm00001d006678\nZm00001d008370\nZm00001d051416\nZm00001d017540\nZm00001d021410")

    
    @dashapp.callback(
        Output('annotation_table', 'children'),
        Input('annotation_gene_list', 'value'))
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
    #                                Symbols
    ###############################################################################

    tab_symbols_content = html.Div([
        dcc.Textarea(
            id='symbols_gene_list',
            value='Textarea content initialized\nwith multiple lines of text',
            style={'width': '100%', 'height': 100},
            ),

        dbc.Button("Show example symbols", color="primary", id="btn_symbols_example", className="mr-1"),

        html.Div(id="symbols_table"),
        ])

    def symbol_select(con, entity_list):
        ls = []
        for entity in entity_list:
            cursorObj = con.cursor()
            cursorObj.execute('''SELECT gene_id, gene_symbol 
                                FROM gene_symbols 
                                WHERE gene_id =  ?  ''', (entity,))
            # (name,) - need the comma to treat it as a single item and not list of letters
            selected = cursorObj.fetchall()
            if selected == []:
                ls.append("")
            else:
                ls.append(selected[0][1])    
        return(ls)

    @dashapp.callback(
        Output('symbols_gene_list', 'value'),
        Input('btn_symbols_example', 'n_clicks'),
        prevent_initial_call=True,)
    def symbols_example(value):
        return("Zm00001d021929\nZm00001d006678\nZm00001d008370\nZm00001d051416\nZm00001d017540\nZm00001d021410")

    @dashapp.callback(
        Output('symbols_table', 'children'),
        Input('symbols_gene_list', 'value'))
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
            output_df = pd.DataFrame({"GeneID": gene_list, "annotation": symbol_select(con, gene_list) })
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
            options=distinct_db_vals(sqnce_path, "gene_families", "genome_id", 
                                    [["All", "all"]]),
            value="all",
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
        if  "all" in genotypes or genotypes == "all":
            genotypes = distinct_db_vals(sqnce_path, "gene_families", "genome_id", return_ls=True) 
        genotypes = str("','".join(genotypes))
        families = str("','".join(families))
        df = pd.read_sql_query("""SELECT protein_id, genome_id, family_name 
                                FROM gene_families
                                WHERE genome_id IN ('{0}') AND family_name IN ('{1}')""".format(genotypes, families), con)
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
        dbc.Button("Show example families", color="primary", id="btn_family_example", className="mr-1"),
        html.Div(id="family_geneIDs_table"),
        ])
    
    def family_gene_select(con, gene_list):
        # Use an input list of genes to find their family assignments
        genotype_ls = []
        family_ls = []
        for gene in gene_list:
            cursorObj = con.cursor()
            cursorObj.execute('''SELECT genome_id, family_name 
                                FROM gene_families
                                WHERE protein_id =  ? ''', (gene,))
            selected = cursorObj.fetchall()
            if selected == []:
                genotype_ls.append("Gene not found")
                family_ls.append("Gene not found")
            else:
                genotype_ls.append(selected[0][0])
                family_ls.append(selected[0][1])
        return([genotype_ls, family_ls])
    
    @dashapp.callback(
        Output('family_gene_list', 'value'),
        Input('btn_family_example', 'n_clicks'),
        prevent_initial_call=True,)
    def family_example(value):
        return("Zm00001d021929\nZm00001d006678\nZm00001d008370\nZm00001d051416\nZm00001d017540\nZm00001d021410")


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
            db_output = family_gene_select(con, gene_list)
            output_df = pd.DataFrame({"GeneID": gene_list, "Genotype ID": db_output[0], "FamilyID": db_output[1] })
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
    #                                Orthogroups
    ###############################################################################
    tab_orthogroups_content = html.Div([
        html.Div(
            className="row", children=[
                html.Div(children=[html.P("Add a gene ID below to get associated orthogroup:", style={'marginTop': 10, 'marginBottom': 10})], 
                    className='six columns', style=dict(width='25%')), 
            ], style=dict(display='flex')),
        
        html.Div(
            className="row", children=[
                html.Div(children=[dcc.Textarea(id='orthogroups_gene_list',
                    value='Add a single gene ID here\n(If multiple gene IDs are added only the first will count)',
                    style={'width': '80%', 'height': 60, 'marginBottom': 10, 'marginTop': 0},
                    )], className='six columns', style=dict(width='40%')),
            ], style=dict(display='flex')),

        html.Div(
            className="row", children=[
                html.Div(children=[dbc.Button("Find orthogroup", 
                                    color="primary", id="btn_submit_orthogroup", className="mr-1",
                                    style={'marginTop': 0, 'marginBottom': 20})], 
                    className='six columns', style=dict(width='25%')), 
            ], style=dict(display='flex')),
                
        html.Div(id="orthogroup_species_dropdown"),
        #html.Div(id="orthogroup_species_copy_genes"),
        html.Div(id="orthogroups_seq_table"),
        ], style=dict(marginLeft=20, marginRight=20))

    @dashapp.callback(
        Output('orthogroup_species_dropdown', 'children'),
        Output('orthogroups_seq_table', 'children'),
        Input('btn_submit_orthogroup', 'n_clicks'),
        State('orthogroups_gene_list', 'value'),
        prevent_initial_call=True,)
    def get_protein_list(_, gene_id):
        con = sqlite3.connect(sqnce_path) # deploy with this
        cursorObj = con.cursor()
        try:
            if "\n" in gene_id:
                gene_id = gene_id.split("\n")[0]
            print("test", gene_id)
            df = pd.read_sql_query("""SELECT orthogroup, genome_id, gene_id
                    FROM orthogroups
                    WHERE orthogroup in (
                            SELECT orthogroup
                            FROM orthogroups
                            WHERE gene_id = '{0}' )""".format(gene_id), con)
            selected = cursorObj.fetchall()
            #if selected == []:
            #    return([html.P("Did not find an orthogroup for the gene."), html.P("")])
        except:
            return([html.P("Something did not work reading the gene list"), html.P("")])
       
        try:
            dropdown_dict = df.groupby("genome_id").count()["orthogroup"].to_dict()
            dropdown_ls = []
            dropdown_ls.append({'label': 'All', 'value': 'all'})
            dropdown_ls.append({'label': 'ZmNAM', 'value': 'ZmNAM'})
            for key in dropdown_dict:
                dropdown_ls.append({'label': key + " ("+ str(dropdown_dict[key]) + ")", 'value': key})

            orthogroup_dropdown_div = html.Div(
                # https://css-tricks.com/wp-content/uploads/2022/02/css-flexbox-poster.png
                className="row", children=[
                    html.Div(children=[dcc.Dropdown(
                        id='orthogroup_select_species_dropdown',
                        options=dropdown_ls,
                        value="all",
                        multi=True,
                        searchable=True
                    ),], className='six columns', style=dict(width='60%')), 
                    html.Div(children=[html.P("Copy selected gene IDs:")],
                        className='six columns', style=dict(width='12%')),
                    html.Div(children=[dcc.Clipboard(id="orthogroups_selected_copy", style={"fontSize":20})], 
                        className='six columns', style=dict(width='5%')),
                    html.Div(children=[dbc.Button("Download Orthogroup Table", 
                            color="primary", id="btn_download_orthogroup", className="mr-1"),
                            dcc.Download(id="download_orthogroup"),],
                            className='six columns', style=dict(width='15%')), 
                ], style=dict(display='flex', justifyContent='flex-start', marginbottom="20px")),

            #orthogroups_table_copy = dcc.Clipboard(id="orthogroups_table_copy", style={"fontSize":20}),

            orthogroup_table = dash_table.DataTable(
                id="orthogroups_table_state",
                columns=[{"name": i, "id": i} for i in df.columns],
                data=df.to_dict('records'),
                style_cell={'textAlign': 'left',
                            'overflow': 'hidden',
                            'maxWidth': 0,},
                editable=True,
                row_deletable=True,)

            return(orthogroup_dropdown_div, orthogroup_table)
            #return([orthogroup_dropdown, orthogroups_table_copy, ])
        except:
            return([html.P("Something did not work reading the gene list"), html.P("")])

    @dashapp.callback(
        Output("orthogroups_selected_copy", "content"),
        Input("orthogroups_selected_copy", "n_clicks"),
        State("orthogroup_select_species_dropdown", "value"),
        State("orthogroups_table_state", "data"),
        prevent_initial_call=True,
    )
    def orthogroup_genes_copy(_, species, df):
        df = pd.DataFrame.from_dict(df)
        gene_groups_ZmNAM = ["ZmB73v4","ZmB97","ZmCML52","ZmCML69","ZmCML103","ZmCML228","ZmCML247","ZmCML277","ZmCML322","ZmCML333","ZmHP301","ZmIl14H","ZmKi3","ZmKi11","ZmKy21","ZmM37W","ZmM162W","ZmMo18W","ZmMs71","ZmNC350","ZmNC358","ZmOh7B","ZmOh43","ZmP39","ZmTx303","ZmTzi8"]
        # If only one item is selected it returns a string instead of a list
        if type(species) == str:
            species = [species]
        
        if "all" in species:
            gene_list = df["gene_id"].to_list()
        else:
            gene_list = df[df["genome_id"].isin(species)]["gene_id"].to_list()

        if "ZmNAM" in species:
            group_gene_list = df[df["genome_id"].isin(gene_groups_ZmNAM)]["gene_id"].to_list()
            gene_list = list(set(gene_list+group_gene_list))
        return ("\n".join(gene_list))

    @dashapp.callback(
        Output("download_orthogroup", "data"),
        Input("btn_download_orthogroup", "n_clicks"),
        State("orthogroups_table_state", "data"),
        prevent_initial_call=True,
    )
    def orthogroup_table_download(_, df):
        df = pd.DataFrame.from_dict(df)
        output = df.to_csv(index=False, header=True, sep="\t", quoting=csv.QUOTE_NONE)
        return dict(content=output, filename="orthogroup.tsv")


    ###############################################################################
    #                                BBHs
    ###############################################################################
    BBHs_query_genotype_list = pd.read_csv(os.path.join(cwd,"init/BBHs_combs.tsv"), sep="\t")
    BBHs_query_genotype_list = list(BBHs_query_genotype_list["query"].unique())
    tab_BBHs_content = html.Div([
        dcc.Dropdown(
            id='dropdown_BBHs_column_select',
            options=[{'label': "Search in Reference Gene IDs (B73v4 and Col-0)", 'value': 'subject_id'},
                     {'label': "Search in All Queried Genomes", 'value': 'query_id'}],
            value="subject_id"
        ),
        dcc.Dropdown(
            id='dropdown_BBHs_query_genotype_select',
            options=[{'label': qg, 'value': qg} for qg in BBHs_query_genotype_list],
            value="All",
            multi=True,
            searchable=True,
        ),
        dcc.Dropdown(
            id='dropbox_show_best_only',
            options=[{'label': "Show Best Blast Hit Only", 'value': 'only'},
                     {'label': "Show Top Best Blast Hits", 'value': 'all'}],
            value="only"
        ),

        dcc.Textarea(
            id='BBHs_gene_list',
            value='Inert gene IDs to search for, one row per gene.',
            style={'width': '100%', 'height': 100},
        ),
        dbc.Button("Show example BBHs", color="primary", id="btn_BBHs_example", className="mr-1"),
        html.Div(id="BBHs_missing"),
        html.Div(id="BBHs_table"),
        ])

    @dashapp.callback(
        Output('BBHs_gene_list', 'value'),
        Output('dropdown_BBHs_column_select', 'value'),
        Output('dropdown_BBHs_query_genotype_select', 'value'),
        Output('dropbox_show_best_only', 'value'),
        Input('btn_BBHs_example', 'n_clicks'),
        prevent_initial_call=True,)
    def BBHs_example(value):
        print("stestfd")
        return "Zm00001d021929", "subject_id", 'All',  "only"


        # selected: variable to select either subject_id or query_id from BBHs table
    @dashapp.callback(
        Output("BBHs_table", "children"),
        Output("BBHs_missing", "children"),
        Input("dropdown_BBHs_column_select", "value"),
        Input("dropdown_BBHs_query_genotype_select", "value"),
        Input("dropbox_show_best_only", "value"),
        Input("BBHs_gene_list", "value"),
        
    )
    def BBHs_select(selected, genotypes, show_best_only, entity_list):
        print(entity_list)
        gene_list = entity_list.split("\n")
        if len(gene_list) > 500:
            gene_list = gene_list[:500]
        if gene_list[-1]=="":
            gene_list = gene_list[:-1]
        
        con = sqlite3.connect(sqnce_path) # deploy with this
        cursorObj = con.cursor()
        gene_list_str = str("','".join(gene_list))
        df = pd.read_sql_query("""SELECT * 
                        FROM BBHs
                        WHERE {0} IN ('{1}')""".format(selected, gene_list_str), con)
        try:
            if "All" not in genotypes:
                df = df[df["query_genome"].isin(genotypes)]
            if show_best_only == "only":
                # This is a stupid function but it seems to work correctly.
                df = df.sort_values(['bit_score'], ascending=False).groupby(["subject_id", "query_genome"]).agg({"bit_score": "first", "query_id": "first",}).reset_index()
            missing_str = ''
            for gene in gene_list:
                if gene not in list(df["subject_id"]):
                    missing_str = missing_str + gene + ", "
            return dash_table.DataTable(
                    columns=[{"name": i, "id": i} for i in df.columns],
                    data=df.to_dict('records'),
                    style_cell={'textAlign': 'left',
                                'overflow': 'hidden',
                                'maxWidth': 0,},
                    editable=True,
                    row_deletable=True,), html.P("The following were not found: "+missing_str),
        except:
            return(html.P("Something did not work returning"),
                   html.P("the Best Blast Hits."))

    ###############################################################################
    #                                Protein Sequences
    ###############################################################################

    tab_protein_content = html.Div([
        dcc.Textarea(
            id='protein_gene_list',
            value='Textarea content initialized\nwith multiple lines of text',
            style={'width': '100%', 'height': 100},
            ),
        dbc.Row(dbc.Col(
            [
            dbc.Button("Download fasta with gene IDs", color="primary", id="btn_download_fasta_geneIDs", className="mr-1"),
            dcc.Download(id="download_fasta_geneIDs"),
            dbc.Button("Download fasta with symbols", color="primary", className="mr-1"),
            dbc.Button("Show example sequences", color="primary", id="btn_protein_example", className="mr-1"),
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
            selected = cursorObj.fetchall()
            if selected == []:
                continue # Skip if cannot find gene
            selected = selected[0] # Otherwise, get the value of the first returned row
            od[selected[0]] = zlib.decompress(selected[1]).decode('utf-8')[:-1] 
        return(od)

    @dashapp.callback(
        Output('protein_gene_list', 'value'),
        Input('btn_protein_example', 'n_clicks'),
        prevent_initial_call=True,)
    def proteins_example(value):
        return("Zm00001d021929\nZm00001d006678\nZm00001d008370\nZm00001d051416\nZm00001d017540\nZm00001d021410")

    @dashapp.callback(
        Output('protein_seq_table', 'children'),
        Input('protein_gene_list', 'value'))
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
            id='promoters_kind_dropdown',
            options=[{'label': "TSS", 'value': "TSS"}, {'label': "ATG", 'value': "ATG"}],
            value="TSS"
            ),
        dcc.Textarea(
            id='promoter_gene_list',
            value='Insert list of genes to\nget a list of promoter sequences',
            style={'width': '100%', 'height': 100},
            ),
        dbc.Row(dbc.Col(
            [
            dbc.Button("Download fasta with gene IDs", color="primary", id="btn_promoter_download_fasta_geneIDs", className="mr-1"),
            dcc.Download(id="download_promoter_fasta_geneIDs"),
            dbc.Button("Download fasta with symbols", color="primary", className="mr-1"),
            dbc.Button("Show example promoters", color="primary", id="btn_promoters_example", className="mr-1"),
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
        Output('promoters_kind_dropdown', 'value'),
        Output('promoter_gene_list', 'value'),
        Input('btn_promoters_example', 'n_clicks'),
        prevent_initial_call=True,)
    def promoter_example(value):
        return "ATG", "Zm00001d021929\nZm00001d006678\nZm00001d008370\nZm00001d051416\nZm00001d017540\nZm00001d021410"


    @dashapp.callback(
        Output('promoter_seq_table', 'children'),
        Input('promoter_gene_list', 'value'),
        Input('promoters_kind_dropdown', 'value'))
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

    ###############################################################################
    #                                Transcriptomics
    ###############################################################################

    tab_transcriptomics_content = html.Div([
        dcc.Textarea(
            id='transcriptomics_gene_list',
            value='Textarea content initialized\nwith multiple lines of text',
            style={'width': '100%', 'height': 100},
            ),
        dbc.Button("Show example transcriptomics", color="primary", id="btn_transcriptomics_example", className="mr-1"),
        html.Div(id="transcriptomics_annotation_table"),
        ])

    def transcriptomic_select(entity_list):
        ls = []
        con = duckdb.connect(database=os.path.join(cwd, "SQNce.duckdb"), read_only=True)
        # Read csv file as pandas dataframe
        query_list = str("','".join(entity_list))
        df = con.execute("""SELECT * FROM PanNAM
                            WHERE target_id IN ('{0}');""".format(query_list)).fetchdf()
        return(df)

    @dashapp.callback(
        Output('transcriptomics_gene_list', 'value'),
        Input('btn_transcriptomics_example', 'n_clicks'),
        prevent_initial_call=True,)
    def symbols_example(value):
        return("Zm00001d021929\nZm00001d006678\nZm00001d008370\nZm00001d051416\nZm00001d017540\nZm00001d021410")


    @dashapp.callback(
        Output('transcriptomics_annotation_table', 'children'),
        Input('transcriptomics_gene_list', 'value'))
    def get_gene_list(value):
        try:
            gene_list = value.split("\n")
            if len(gene_list) > 500:
                gene_list = gene_list[:500]
            if gene_list[-1]=="":
                gene_list = gene_list[:-1]
        except:
            return(html.P("Something did not work reading the gene list"))
        try:
            output_df = transcriptomic_select(gene_list)
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
        