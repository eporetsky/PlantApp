from numpy import Inf


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
    import subprocess
    import zlib
    import base64
    import io

    from collections import OrderedDict

    import dash
    from dash import no_update # https://community.plotly.com/t/error-expected-the-output-type-to-be-a-list-or-tuple-but-got-none/34744/6
    from flask import Flask, send_from_directory
    from urllib.parse import quote as urlquote

    from io import BytesIO
    from PIL import Image

    from Bio import SeqIO
    from Bio import Phylo
    from Bio import AlignIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from PIL import Image


    # import numpy as np
    from io import BytesIO
    import base64

    from flask_login import current_user
    from flask import session

    # import pandas as pd
    from plotly.tools import mpl_to_plotly
    import plotly.express as px
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    # Set the cwd and SQNce.db path to be platform independent  
    if "plantapp" not in os.getcwd().lower():
        cwd = "/home/eporetsky/plantapp" # for server hosting
    else:
        cwd = os.getcwd() # for personal computer
    sqnce_path = os.path.join(cwd, "SQNce.db")
    
    @dashapp.callback(Output('tab_content', 'children'),
                      Input('url', 'pathname'))
    def display_page_content(pathname):
        #print(session.get('username', None))
        pathname = pathname.split("/")[-1]
        if pathname == "simple_tree":
            return tab_simple_tree
        if pathname == "genome_graph":
            return tab_genome_graph
        if pathname == "mapping_summary":
            return tab_mapping_summary



    ###############################################################################
    #                                Simple Tree
    ###############################################################################


    tab_simple_tree = html.Div([
        dcc.Textarea(
                id='simple_tree_gene_list',
                value='Paste gene list',
                style={'width': '100%', 'height': 100},
                ),
        html.Button('Prepare alignment and tree', id='simple_tree_aln_button', type='submit'),
        html.Div(id='simple_tree_tree_figure'),
        #html.Div([html.Img(id = 'simple_tree_tree_figure', src = '')], id='plot_div'),
        ])

    def simple_tree_write_fasta(con, entity_list, session_id):
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
            record = SeqRecord(Seq(zlib.decompress(selected[1]).decode(encoding='UTF-8')[:-1]), id=selected[0], name="", description="")
            od[selected[0]] = record
        with open(os.path.join(cwd,"pages", "Apps", "download", f"{session_id}_selected.fasta"), 'w') as handle:
            SeqIO.write(od.values(), handle, 'fasta')

    @dashapp.callback(
        Output("simple_tree_tree_figure", "children"),
        Input('simple_tree_aln_button', 'n_clicks'),
        State('simple_tree_gene_list', 'value'))
    def alignment_update(clicks, value):
        import random
        import string
        session_id = ''.join(random.SystemRandom().choice(string.ascii_letters + string.digits) for _ in range(10))

        # Not sure why this callback is triggered when loading the app. This is a way out.
        if callback_context.triggered[0]["prop_id"] == ".":
            return(html.P("Insert list of genes to generate the fasta file."))

        gene_list = value.split("\n")
        if len(gene_list) > 50:
            gene_list = gene_list[:50]
        if gene_list[-1]=="":
            gene_list = gene_list[:-1]
        try:
            con = sqlite3.connect(sqnce_path)
            simple_tree_write_fasta(con, gene_list, session_id)
            #print("ran protein")
            con.close()
            # uses a separate function to generate the download link
        except:
            return(html.P("Something did not work writing the fasta file"))

        ctx = dash.callback_context
        if (not ctx.triggered and not ctx.triggered[0]['value'] == 0):
            return (no_update)
        if clicks is not None:
            try:
                app_path = os.path.join(cwd,"pages", "Apps")
                print("running famsa") # Needs: chmod a+x famsa
                subprocess.run([app_path+"/./famsa "+app_path+f"/download/{session_id}_selected.fasta "+app_path+f"/download/{session_id}_selected.aln", "arguments"], shell=True)
                print("running clipkit") # Needs pip or conda install, make sure I didn't use binary when deployed
                subprocess.run(["clipkit "+app_path+f"/download/{session_id}_selected.aln", "arguments"], shell=True)
                print("running fasttree") # Needs: chmod a+x fasttree
                subprocess.run([app_path+"/./fasttree  "+app_path+f"/download/{session_id}_selected.aln.clipkit > "+app_path+f"/download/{session_id}_selected.tree", "arguments"], shell=True)
                ########################################
                # Now we should have a tree and an alignment file in the download folder
                # In 3.x, there is no longer a .next method; use the built-in function:
                tree = next(Phylo.parse(app_path+f'/download/{session_id}_selected.tree', 'newick'))
                tree.root_at_midpoint() # operates in-place
                align = AlignIO.read(app_path+f"/download/{session_id}_selected.aln", "fasta")
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, len(align)), constrained_layout=True)
                #fig.suptitle('Horizontally stacked subplots')
                # https://stackoverflow.com/questions/29419973/in-python-how-can-i-change-the-font-size-of-leaf-nodes-when-generating-phylogen
                ########################################
                line_width = 0.1 # sometimes if the line is not thick enought it might not appear in the %inline
                box_width = 0.5
                y_pad = 0.1
                #add rectangle to plot
                def names(terminal):
                    return terminal.name
                phylo_ordered_list = list(map(names, tree.get_terminals()))

                for aln in align:
                    aln_count = len(align) - phylo_ordered_list.index(aln.id) - 1
                    ax2.add_patch(Rectangle((0, aln_count - line_width/2), len(align[0].seq), line_width))
                    rect_ls = alignment_to_range_list(aln.seq)
                    # Read every other value in the rect_ls> Paired values correspond to start and end
                    for rect in range(0, len(rect_ls), 2):
                        ax2.add_patch(Rectangle((rect_ls[rect], aln_count-box_width/2), rect_ls[rect+1] - rect_ls[rect], box_width))
                    aln_count -= 1
                ax2.set(yticklabels=[])
                ax2.set(ylabel=None)
                ax2.tick_params(left=False)  # remove the ticks
                ax2.margins(y=(2/(len(align)*3))) # I don't know why it works but close enough
                ax2.plot(0)
                
                ########################################
                # Keep this reverse order because I Phylo.draw uses plt.show() (I assume) which prevents ax2 from showing up
                matplotlib.rc('font', size=8)
                #plt.figure(figsize=(20, 6))
                Phylo.draw(tree, axes=ax1, do_show=False)
                #ax.add_patch(Rectangle((box_width, box_width), 2, 6))
                fig.savefig(app_path+f"/download/{session_id}_selected_tree_aln.png")
                ########################################
                return(html.Img(src = fig_to_uri(fig)))
            except:
                return(html.P("Something did not work creating the image."))

            

    def alignment_to_range_list(alignment):
        # function takes a single sequence from an alignment
        # returns a 1-d list of paired start and end coordinates for each aligned segment
        count = -1 # because I use continue I moved the counter to the beginning of the loop
        streak = False # keeps track of streaks of non-gapped characters in the loop progression
        range_list = []
        for s in alignment:
            count += 1
            if s == "-":
                if streak == True:
                    range_list.append(count)
                    streak = False
                continue
            elif streak == True:
                continue
            else:
                range_list.append(count)
                streak = True
        # When alignment doesn't end with a gap finish the list
        if len(range_list)%2:
            range_list.append(len(alignment))
        return range_list

    


    ###############################################################################
    #                                Genome Graph
    ###############################################################################

    tab_genome_graph = html.Div([
        # Get a list of unique genotypes from the db
        dcc.Dropdown(
            id='genome_graph_genotypes_dropdown',
            options=[
                {'label': 'B73v4', 'value': "B73v4"},
                {'label': 'Arabidopsis', 'value': "Arabidopsis"},
            ],
            value="None"
        ),

        dcc.Textarea(
                id='genome_graph_list',
                value='Paste gene list or coordiantes (tab separated Chromsome and Location',
                style={'width': '100%', 'height': 100},
                ),
        html.Button('Prepare genome graph figure', id='genome_graph_button', type='submit'),
        html.Div(id='genome_graph_figure'),
        #html.Div([html.Img(id = 'simple_tree_tree_figure', src = '')], id='plot_div'),
    ])

    # Query to find neighboring genes
    def get_gene_coordinates(genotype, gene_id):
        con = sqlite3.connect(sqnce_path) # deploy with this
        cursorObj = con.cursor()
        df = pd.read_sql_query('''SELECT * 
                        FROM gene_coordinates 
                        WHERE genotype_id = "{0}"
                        AND gene_id = "{1}"
                        '''.format(genotype, gene_id), con)
        return(df)

    @dashapp.callback(
        Output("genome_graph_figure", "children"),
        Input('genome_graph_button', 'n_clicks'),
        Input("genome_graph_genotypes_dropdown", 'value'),   
        State('genome_graph_list', 'value'))
    def genome_graph_update(clicks, genotype, gene_list):
        import random
        import string
        session_id = ''.join(random.SystemRandom().choice(string.ascii_letters + string.digits) for _ in range(10))
        app_path = os.path.join(cwd,"pages", "Apps")
        # Not sure why this callback is triggered when loading the app. This is a way out.
        if callback_context.triggered[0]["prop_id"] == ".":
            return(html.P("Insert list of genes to generate the fasta file."))

        gene_list = gene_list.split("\n")
        if len(gene_list) > 500:
            gene_list = gene_list[:500]
        if gene_list[-1]=="":
            gene_list = gene_list[:-1]
        try:
            # Get the gene coordinates for the genome graph
            coords = pd.DataFrame()
            for gene in gene_list:
                coords = pd.concat([coords, get_gene_coordinates(genotype, gene)])
        except:
            return(html.P("Something did not work writing the fasta file"))
        ctx = dash.callback_context
        if (not ctx.triggered and not ctx.triggered[0]['value'] == 0):
            return (no_update)
        if clicks is not None:
            try:
                from reportlab.lib.units import cm
                from reportlab.graphics import renderPM
                from Bio import SeqIO
                from Bio.Graphics import BasicChromosome
                df = pd.read_csv(os.path.join(app_path, "chromosome.tsv"), sep="\t")
                selected = df[df["genotype_id"]==genotype]
                entries = selected[["chr", "len"]].values.tolist()
                max_len = selected["len"].max()
                telomere_length = 0
                chr_diagram = BasicChromosome.Organism()
                chr_diagram.page_size = (29.7 * cm, 21 * cm)  # A4 landscape

                for index, (name, length) in enumerate(entries):
                    
                    # For each chrosome generate a feature list of tuples (start, start, strand, name, color,)
                    features = []
                    for row in coords[coords["gene_chr"]==name].values:
                        features.append((row[3], row[3], None, row[0], "black",))
                    cur_chromosome = BasicChromosome.Chromosome(name)
                    # Set the scale to the MAXIMUM length plus the two telomeres in bp,
                    # want the same scale used on all five chromosomes so they can be
                    # compared to each other
                    cur_chromosome.scale_num = max_len + 2

                    # Record an Artemis style integer color in the feature's qualifiers,
                    # 1 = Black, 2 = Red, 3 = Green, 4 = blue, 5 =cyan, 6 = purple
                    
                    # The features can either be SeqFeature objects, or tuples of values: start (int), 
                    # end (int), strand (+1, -1, O or None), label (string), ReportLab color (string or object), 
                    # and optional ReportLab fill color.
                    
                    # This will it to every chromsome. Make chromsome specific lists
                    
                    # Add an opening telomere
                    start = BasicChromosome.TelomereSegment()
                    start.scale = telomere_length
                    cur_chromosome.add(start)
                    
                    # Add a body - again using bp as the scale length here.
                    body = BasicChromosome.AnnotatedChromosomeSegment(length, features)
                    body.scale = length
                    cur_chromosome.add(body)
                    
                    # Add a closing telomere
                    end = BasicChromosome.TelomereSegment(inverted=True)
                    end.scale = telomere_length
                    cur_chromosome.add(end)
                    
                    # This chromosome is done
                    chr_diagram.add(cur_chromosome)
                pdf_path = os.path.join(app_path, "download", f"{session_id}_genome_graph.pdf")
                chr_diagram.draw(pdf_path, "Genome Diagram")
                fig = chr_diagram.draw(None, "Genome Diagram")
                #print(os.path.join(app_path,"download", f"{session_id}_genome_graph.png"))

                #print(os.path.join(cwd,"download", "Apps", f"{session_id}_genome_graph.png"))
                # This saves it as a png file, which I think that I don't need to draw it on page
                #fig = renderPM.drawToFile(fig, pdf_path.replacte("pdf", "png"), 'PNG')
                print("final test")

                with open(pdf_path, 'rb') as pdf:
                    pdf_data = base64.b64encode(pdf.read()).decode('utf-8')


                #return(html.Iframe(id="embedded-pdf", src=pdf_path))
                return(
                    html.ObjectEl(
                        # To my recollection you need to put your static files in the 'assets' folder
                        data='data:application/pdf;base64,'+ pdf_data,
                        type="application/pdf",
                        style={"width": "100%", "height": "1000px"}
                         ),
                )
                    
                    
                #return(html.Img(src = fig_to_uri(fig)))
            except:
                return(html.P("Something did not work creating the image."))


    def fig_to_uri(in_fig, close_all=True, **save_args):
        """
        Save a figure as a URI
        :param in_fig:
        :return:
        """
        out_img = BytesIO()
        in_fig.savefig(out_img, format='png', **save_args)
        if close_all:
            in_fig.clf()
            plt.close('all')
        out_img.seek(0)  # rewind file
        encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
        return "data:image/png;base64,{}".format(encoded)

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

    def return_column_if(db, table, return_column, check_column, condition, custom_vals=[]):
        ls = [{ 'label': label, 'value': val} for label, val in custom_vals]
        con = sqlite3.connect(db) # deploy with this
        cursorObj = con.cursor()
        distinct_df = pd.read_sql_query('''SELECT {0} 
                                        FROM {1}
                                        WHERE {2}=="{3}"'''.format(return_column, table, check_column, condition), con)
        for name in distinct_df[return_column]:
            ls.append({'label': name, 'value': name})
        return(ls)

    def parse_contents(contents, filename, date):
        if filename.split(".")[-1] in ["png", "jpg", "jpeg"]:
            return html.Div([
                html.H5(filename),
                # HTML images accept base64 encoded strings in the same format
                # that is supplied by the upload
                
                html.Img(src=contents, style={'width': '80%', 'height': 'auto'}),
                html.Hr(),
                html.Div('Raw Content'),
                html.Pre(contents[0:200] + '...', style={
                    'whiteSpace': 'pre-wrap',
                    'wordBreak': 'break-all'
                })
                ], style = {'max-width': '1200px'})

        if filename.split(".")[-1] == "csv":
            contents = contents.split(',')[-1]
            decoded = base64.b64decode(contents)
            # Not sure if there's a better way to remove the file header "data:application/vnd.ms-excel;base64"
            #print(decoded[0])
            #decoded = decoded[0].split("base64,")[-1]
            df = pd.read_csv(io.StringIO(decoded.decode('utf-8')))
            df = df.sort_values(df.columns[-1])
            return html.Div([
                html.H5(filename),
                dash_table.DataTable(
                    id='datatable-interactivity',
                    columns=[{"name": i, "id": i, "deletable": True, "selectable": True} for i in df.columns],
                    data=df.to_dict('records'),
                    editable=False,
                    filter_action="native",
                    sort_action="native",
                    sort_mode="multi",
                    column_selectable="single",
                    row_selectable="multi",
                    row_deletable=True,
                    selected_columns=[],
                    selected_rows=[],
                    page_action="native",
                    page_current= 0,
                    page_size= 10,
                ),

            ], style = {'max-width': '1200px'})

    ###############################################################################
    #                                Mapping Summary
    ###############################################################################
    
    tab_mapping_summary = html.Div([
        html.Br(),
        html.Div([
            dbc.Row([
                dbc.RadioItems(
                    id="mapping_summary_select_menu",
                    className="btn-group",
                    inputClassName="btn-check",
                    labelClassName="btn btn-outline-primary",
                    labelCheckedClassName="active",
                    options=[
                        {"label": "Upload Data", "value": 1},
                        {"label": "Manage Data", "value": 2},
                        {"label": "Annotate Results", "value": 3},
                        {"label": "Candidates Table", "value": 4},
                        {"label": "Candidates Info", "value": 5},
                    ], value=5,
                    ),
                ], align="end",
                ),
        ], style={"height": "50px"},
        ),
        
        html.Div(id="mapping_summary_select_menu_output")
    ])

    @dashapp.callback(Output("mapping_summary_select_menu_output", "children"), 
                      Input("mapping_summary_select_menu", "value"))
    def mapping_summary_display_page(value):
        #print(value)
        if value == 1:
            return(tab_mapping_summary_upload)
        if session.get('username', None)==None:
            return(tab_mapping_login_request)
        if value == 2:
            return(tab_mapping_summary_manage)
        if value == 3:
            return(tab_mapping_summary_annotate)
        #print(session.get('username', None))
        if value == 4:
            return(tab_mapping_summary_table)
        if value == 5:
            return(tab_mapping_summary_info)

    tab_mapping_summary_upload = html.Div([
        html.P("This is the upload tab."),
        dcc.Upload(
            id='upload-image',
            children=html.Div([
                'Drag and Drop or ',
                html.A('Select Files')
            ]),
            style={
                'width': '100%',
                'height': '60px',
                'lineHeight': '60px',
                'borderWidth': '1px',
                'borderStyle': 'dashed',
                'borderRadius': '5px',
                'textAlign': 'center',
                'margin': '10px'
            },
            # Allow multiple files to be uploaded
            multiple=True
        ),
        html.Div(id='output_image_upload'),
        html.Div(id='div_selected_table', style = {'max-width': '1200px'}),
    ])

    @dashapp.callback(
              Output('output_image_upload', 'children'),
             #Output('div_selected_table', 'children'),
              Input('upload-image', 'contents'),
              State('upload-image', 'filename'),
              State('upload-image', 'last_modified'))
    def update_output(list_of_contents, list_of_names, list_of_dates):
        if list_of_contents is not None:
            children = [
                parse_contents(c, n, d) for c, n, d in
                zip(list_of_contents, list_of_names, list_of_dates)]
            return children

    @dashapp.callback(
        Output('div_selected_table', 'children'),
        Input('datatable-interactivity', 'selected_rows'),
        State('datatable-interactivity', 'data'), 
        prevent_initial_call=True)
    def print_value(selected_rows,rows):
        con = sqlite3.connect(sqnce_path) # deploy with this
        rows = pd.DataFrame.from_dict(rows)
        #print(df)
        #print(selected_rows)
        rows = rows.iloc[selected_rows, :]
        #print(df)
        coordinate_list = rows[["CHROM", "POS"]].values
        # Check the neighboring genes for each row in the input data
        df = pd.DataFrame()
        for row in coordinate_list:
            df = pd.concat([df, get_SNP_neighbors("B73v4", str(row[0]), int(row[1]), 100000)])
        # Since some genes can overlap, remove duplicated values and eep first
        
        df = df.drop_duplicates(subset=["gene_id"]).reset_index()

        # Add a column of gene annotations when available
        df["annotation"] = multiple_row_select(con, "gene_annotations", df["gene_id"].to_list())
        df = df.drop(["index", "genotype_id", "gene_orientation"], axis=1)


        return dash_table.DataTable(
                    id='selected_table',
                    columns=[{"name": i, "id": i, "deletable": False, "selectable": False} for i in df.columns],
                    data=df.to_dict('records'),
                    editable=False,
                    filter_action="native",
                    sort_action="native",
                    sort_mode="multi",
                    column_selectable="single",
                    row_selectable="multi",
                    row_deletable=True,
                    selected_columns=[],
                    selected_rows=[],
                    page_action="native",
                    page_current= 0,
                    page_size= 10,
                )

    ###########################################################################
    ###########################################################################
    tab_mapping_login_request = html.Div([
        html.P("Please login for access.")
    ])

    ###########################################################################
    ###########################################################################
    tab_mapping_summary_manage = html.Div([
        html.P("This tab can be used to see how many experiments are available for the user."),
        # Get a list of unique genotypes from the db
        dcc.Dropdown(
            id='coordinates_genotypes_dropdown',
            options=distinct_db_vals(sqnce_path, "mapping_traits", "trait", [["All", "all"]]),
            value="all"
        ),
        dcc.Dropdown(
            id='coordinates_genotypes_dropdown',
            options=[{'label': "All", 'value': "all"},
                     {'label': "Not annotated", 'value': "unannotated"},
                     {'label': "Ambiguous", 'value': "ambiguous"},
                     {'label': "Annotated", 'value': 'annotated'}],
            value="unannotated"
        ),
    ])

    ###########################################################################
    ###########################################################################

    tab_mapping_summary_annotate = html.Div([
        dcc.Store(id='store_plot_rowid_list'),
        dcc.Store(id='store_current_plot_list_index'),
        html.P("This tab will help access all the mapping candidates."),
        html.Div(id="coordinates_genotypes_dropdown_div",
                 style = {"width": "95%"}),
        # Get a list of unique genotypes from the db
        dcc.Dropdown(
            id='mapping_experiments_dropdown',
            options=distinct_db_vals(sqnce_path, "mapping_traits", "experiment",),
                #[["All","all"], ["Complete","complete"],["Skipped","skipped"],["Unannotated","unannotated"]]),
            value="Select Experiment Name",
            style = {"width": "95%"},
        ),
        html.Div(id='div_mapping_experiments_dropdown', style = {"width": "95%"}), # 'max-width': '95%' 1200px
        html.P(),
        html.Div([dbc.ButtonGroup([
            dbc.Button("Previous plot", id="parser_previous_plot", outline=False),
            dbc.Button("Next plot", id="parser_skip_plot", outline=False),
            dbc.Button("Accept results", id="parser_accept_results", outline=False, color="warning"),
            ], size="lg", style={"width": "95%"}),
        ], style = {"width": "95%"}),
        html.P(),
        html.Div([
            dbc.ButtonGroup([
                dbc.Button("Skip", id="score_skip", outline=True, color="primary"),
                dbc.Button("1", id="score_1", outline=True, color="primary"),
                dbc.Button("2", id="score_2", outline=True, color="primary"),
                dbc.Button("3", id="score_3", outline=True, color="primary"),
                dbc.Button("4", id="score_4", outline=True, color="primary"),
                dbc.Button("5", id="score_5", outline=True, color="primary"),
            ], size="lg", style={"width": "95%"}),
        ], id="mapping_score_dropdown_div", style = {"width": "95%"}),
        html.P(),
        html.Div(id='div_mapping_summary_parser', style = {"width": "95%"}),
        html.P(),
        html.Div(id="mapping_parser_candidates_table", style = {"width": "95%"}),
        html.Div(id='mapping_empty_placeholder', style = {"width": "95%"}),

        # To start I can do a dropboxwith "all" and then groupby and count
        # And otherwise sort by distinct symbols, gene IDs
    ])

    @dashapp.callback(
        Output('div_mapping_experiments_dropdown', 'children'),
        Input('mapping_experiments_dropdown', 'value'),
        prevent_initial_call=True)
    def coordinates_genotypes_dropdown_div(selected_experiments):
        options = distinct_db_vals(sqnce_path, 'mapping_traits', 'processed', [["All", "all"]])
        options += return_column_if(sqnce_path, 'mapping_traits', 'trait', 'experiment', selected_experiments)
        return html.Div([dcc.Dropdown(
            id='select_mapping_traits',
            options=options,
            value="all",
            style = {"width": "95%"},)
        ])

    @dashapp.callback(
        [Output('score_skip', 'outline'), Output('score_skip', 'color'),
         Output('score_1', 'outline'), Output('score_1', 'color'),
         Output('score_2', 'outline'), Output('score_2', 'color'),
         Output('score_3', 'outline'), Output('score_3', 'color'),
         Output('score_4', 'outline'), Output('score_4', 'color'),
         Output('score_5', 'outline'), Output('score_5', 'color')],
        [Input('score_skip', 'n_clicks'),
        Input('score_1', 'n_clicks'),
        Input('score_2', 'n_clicks'),
        Input('score_3', 'n_clicks'),
        Input('score_4', 'n_clicks'),
        Input('score_5', 'n_clicks'),
        Input('store_plot_rowid_list', 'data'),
        Input('store_current_plot_list_index', 'data')],
        prevent_initial_call=True)
    def mapping_score_dropdown_div(skip,s1,s2,s3,s4,s5,rowid_list,current_rowid):
        # I store the score values as 0,1,2,3,4,5 for simplicity
        con = sqlite3.connect(sqnce_path)
        if rowid_list != None:
            rowid = rowid_list[current_rowid]
        cursorObj = con.cursor()
        context = [p["prop_id"] for p in dash.callback_context.triggered][0]
        #print("context:",context, current_rowid)
        default = [True,"primary"]*6
        if context == "score_skip.n_clicks":
            default[0], default[1] = False, "success"
            cursorObj.execute("UPDATE mapping_traits SET score=0 WHERE ROWID={0}".format(rowid))
        if context == "score_1.n_clicks":
            default[2], default[3] = False, "success"
            cursorObj.execute("UPDATE mapping_traits SET score=1 WHERE ROWID={0}".format(rowid))
        if context == "score_2.n_clicks":
            default[4], default[5] = False, "success"
            cursorObj.execute("UPDATE mapping_traits SET score=2 WHERE ROWID={0}".format(rowid))
        if context == "score_3.n_clicks":
            default[6], default[7] = False, "success"
            cursorObj.execute("UPDATE mapping_traits SET score=3 WHERE ROWID={0}".format(rowid))
        if context == "score_4.n_clicks":
            default[8], default[9] = False, "success"
            cursorObj.execute("UPDATE mapping_traits SET score=4 WHERE ROWID={0}".format(rowid))
        if context == "score_5.n_clicks":
            default[10], default[11] = False, "success"
            cursorObj.execute("UPDATE mapping_traits SET score=5 WHERE ROWID={0}".format(rowid))
        if context == "store_current_plot_list_index.data" or context == "store_plot_rowid_list.data":
            df = pd.read_sql_query('''SELECT score 
                FROM mapping_traits
                WHERE rowid={0}'''.format(rowid), con)
            #print(df.values.tolist())
            df_val = df.values.tolist()[0][0]
            #print("val:", df_val)
            if df_val == None: None
            else: default[df_val*2], default[df_val*2+1] = False, "success"
        con.commit()
        con.close()
        return(default)

        
   
    @dashapp.callback(
        [Output('store_plot_rowid_list', 'data'),
        Output('store_current_plot_list_index', 'data')],
        [Input('select_mapping_traits', 'value'),
        Input('parser_previous_plot', 'n_clicks'),
        Input('parser_skip_plot', 'n_clicks'),
        Input('parser_accept_results', 'n_clicks'),],
        [State('mapping_experiments_dropdown', 'value'),
        State('store_plot_rowid_list', 'data'),
        State('store_current_plot_list_index', 'data',)],
        prevent_initial_call=True)
    def store_plot_rowid_list(selected_traits, previous, skip, submit, 
                              mapping_experiment ,rowid_list, current_index):
        context = [p["prop_id"] for p in dash.callback_context.triggered][0]
        if rowid_list != None:
            rowid = rowid_list[current_index]

        con = sqlite3.connect(sqnce_path)
        
        if context == "parser_previous_plot.n_clicks":
            if current_index == 0:
                return(dash.no_update, dash.no_update)
            else:
                return(dash.no_update, current_index-1)
        if context == "parser_skip_plot.n_clicks":
            cursorObj = con.cursor()
            cursorObj.execute('UPDATE mapping_traits SET processed="Skipped" WHERE ROWID={0}'.format(rowid))
            con.commit()
            if current_index==len(rowid_list)-1:
                return(dash.no_update, dash.no_update)
            else:
                return(dash.no_update, current_index+1)
        if context == "parser_accept_results.n_clicks":
            del rowid_list[current_index]
            return(rowid_list, current_index)

        
        if selected_traits == "all":
            df = pd.read_sql_query('''SELECT rowid 
                FROM mapping_traits
                WHERE experiment="{0}"'''.format(mapping_experiment), 
                con)
        elif selected_traits == "Unannotated":
            df = pd.read_sql_query('''SELECT rowid 
                FROM mapping_traits
                WHERE processed="Unannotated" AND experiment="{0}"'''.format(mapping_experiment),
                con) # IS NULL
        elif selected_traits == "Complete":
            df = pd.read_sql_query('''SELECT rowid 
                FROM mapping_traits
                WHERE processed="Complete" AND experiment="{0}"'''.format(mapping_experiment),
                con)
        elif selected_traits == "Skipped":
            df = pd.read_sql_query('''SELECT rowid 
                FROM mapping_traits
                WHERE processed="Skipped" AND experiment="{0}"'''.format(mapping_experiment),
                con)
        else:
            df = pd.read_sql_query('''SELECT rowid 
                FROM mapping_traits
                WHERE trait="{0}" AND experiment="{1}"'''.format(selected_traits,mapping_experiment), 
                con)
        #print(df["rowid"].values.tolist())
        print(df["rowid"].values.tolist())
        return(df["rowid"].values.tolist(), 0)

    def pil_to_b64(im, enc_format="png", **kwargs):
        """
        Converts a PIL Image into base64 string for HTML displaying
        :param im: PIL Image object
        :param enc_format: The image format for displaying. If saved the image will have that extension.
        :return: base64 encoding
        """
        buff = BytesIO()
        im.save(buff, format=enc_format, **kwargs)
        encoded = base64.b64encode(buff.getvalue()).decode("utf-8")

        return encoded
        

    @dashapp.callback(
        Output('div_mapping_summary_parser', 'children'),
        [Input('store_plot_rowid_list', 'data'),
        Input('store_current_plot_list_index', 'data')],
        prevent_initial_call=True)
    def mapping_parser_div(rowid_list, current_index):
        # https://dash.gallery/dash-image-annotation/
        
        #print(rowid_list,current_index)
        if rowid_list != None:
            rowid = rowid_list[current_index]

        con = sqlite3.connect(sqnce_path)
        cursorObj = con.cursor()

        cursorObj.execute("""SELECT plot, trait
                            FROM mapping_traits
                            WHERE rowid= ? """, (rowid,))
        fig, trait = cursorObj.fetchall()[0]
        
        df = pd.read_sql_query('''SELECT * 
                                  FROM mapping_results
                                  WHERE trait = "{0}"'''.format(trait), con)
        df = df.sort_values(df.columns[-1])
        df = df.drop_duplicates()

        picture_stream = io.BytesIO(fig)
        picture = Image.open(picture_stream)

        image_annotation_card = dbc.Row([
            dbc.Col([html.Img(className="image", src="data:image/png;base64, " + pil_to_b64(picture))],
                width={"size": 4, "offset": 0}),
            dbc.Col([dash_table.DataTable(
                id='mapping_parser_results_table',
                columns=[{"name": i, "id": i, "deletable": False, "selectable": False} for i in df.columns],
                data=df.to_dict('records'),
                editable=False,
                #filter_action="native",
                sort_action="native",
                sort_mode="multi",
                #column_selectable="single",
                row_selectable="multi",
                row_deletable=False,
                #selected_columns=[],
                selected_rows=[0],
                page_action="native",
                page_current= 0,
                page_size= 10,
                hidden_columns=["experiment", "trait", "ref", "alt"],
                css=[{"selector": ".show-hide", "rule": "display: none"}], # removes "Toggle Columns" button

            )], width={"size": 6, "offset": 2, "order": "last"}),
        ]),
        return(image_annotation_card)

    @dashapp.callback(
        Output('mapping_parser_candidates_table', 'children'),
        Input('mapping_parser_results_table', 'selected_rows'),
        State('mapping_parser_results_table', 'data'), 
        prevent_initial_call=True)
    def mapping_parser_candidates_table(selected_rows,rows):
        empty_dash_table = dash_table.DataTable(
                columns=[{"name": i, "id": i, "deletable": False, "selectable": False} for i in ["no neighboring genes found"]],
            )
        
        #print("SELECTED ROWS:", selected_rows)
        if len(selected_rows) == 0:
            #print("TEST TEST TEST")
            return empty_dash_table
        con = sqlite3.connect(sqnce_path)
        rows = pd.DataFrame.from_dict(rows)
        rows = rows.iloc[selected_rows, :]
        coordinate_list = rows[["chrom", "pos", "pval", "trait", "experiment", "effect", "snp"]].values
        df = pd.DataFrame()
        are_all_empty = True
        for row in coordinate_list:
            tmp_df = get_SNP_neighbors("B73v4", str(row[0]), int(row[1]), 100000)
            if tmp_df.empty:
                continue
            are_all_empty = False

            tmp_df.insert(0, "pval",row[2])
            tmp_df.insert(0, "trait",row[3])
            tmp_df.insert(0, "experiment",row[4])
            tmp_df.insert(0, "effect",row[5])
            tmp_df.insert(0, "snp",row[6])
            df = pd.concat([df, tmp_df])
        
        if are_all_empty:
            return empty_dash_table
        df = df.drop_duplicates(subset=["gene_id"]).reset_index()


        # Add a column of gene annotations when available
        df["annotation"] = multiple_row_select(con, "gene_annotations", "gene_annotation", df["gene_id"].to_list())
        df["symbols"] = multiple_row_select(con, "gene_symbols", "gene_symbol", df["gene_id"].to_list())
        df = df.drop(["index", "genotype_id", "gene_orientation"], axis=1)
        
                # Re-order columns so it matches the SQL INSERT function
        df = df[["gene_id","symbols","annotation","gene_start","gene_end","distance",
                  "experiment","trait","snp","chrom","pos","pval","effect"]]
        df = df.sort_values("distance")

        return dash_table.DataTable(
                    id='selected_candidate_table',
                    columns=[{"name": i, "id": i, "deletable": False, "selectable": False} for i in df.columns],
                    data=df.to_dict('records'),
                    editable=False,
                    #filter_action="native",
                    sort_action="native",
                    sort_mode="multi",
                    #column_selectable="single",
                    row_selectable="multi",
                    row_deletable=False,
                    #selected_columns=[],
                    selected_rows=[],
                    page_action="native",
                    page_current= 0,
                    page_size= 10,
                    hidden_columns=["pval", "effect","query", "trait", "experiment","chrom","pos"], #"experiment", "trait",
                    css=[{"selector": ".show-hide", "rule": "display: none"}], # removes "Toggle Columns" button 
                )
        
    @dashapp.callback(
        Output('mapping_empty_placeholder', 'children'),
        Input('parser_accept_results', 'n_clicks'   ),
        [State('selected_candidate_table', 'data'),
        State('selected_candidate_table', 'selected_rows'),
        State('store_plot_rowid_list', 'data'),
        State('store_current_plot_list_index', 'data',)],
        prevent_initial_call=True)
    def mapping_parser_candidates_table(accept_rows, rows, selected_rows, rowid_list, current_index):
        #if len(selected_rows) == 0:
        #    return(html.P("No candidate genes were selected."))
        rows = pd.DataFrame.from_dict(rows)
        rows = rows.iloc[selected_rows, :]
        
        con = sqlite3.connect(sqnce_path)
        rowid = rowid_list[current_index]
        cursorObj = con.cursor()
        # Can't use processed="processed"
        cursorObj.execute('UPDATE mapping_traits SET processed="Complete" WHERE ROWID={0}'.format(rowid))
        cursorObj.execute('UPDATE mapping_traits SET num_candidates={0} WHERE ROWID={1}'.format(len(rows), rowid))
        con.commit()
        #con.close()
        
        gene_list_output = []
        for row in rows.values:
            cursorObj = con.cursor()
            cursorObj.execute('''SELECT *
                FROM mapping_candidates
                WHERE gene_id="{0}" AND trait="{1}"'''.format(row[0], row[7]))
            selected = cursorObj.fetchall()
            # Only add new candidate if gene/trait combination doesn't exist
            if selected==[]:
                gene_list_output.append(row[0])
                cursorObj.execute("""INSERT INTO mapping_candidates(
                         gene_id, 
                         gene_symbol, 
                         gene_annotation,
                         gene_start,
                         gene_end,
                         distance,
                         experiment,
                         trait,
                         snp,
                         chrom,
                         pos,
                         pval,
                         effect) 
                         VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)""", row)
                con.commit()
        
        return html.P("adding the following genes: "+", ".join(gene_list_output))
    ###########################################################################
    ###########################################################################
    
    tab_mapping_summary_table = html.Div([
        dbc.Label("This tab will help access all the mapping candidates."),
        # To start I can do a dropboxwith "all" and then groupby and count
        # And otherwise sort by distinct symbols, gene IDs
        dbc.Row([
            dbc.Col(html.P("Select which experiments to show:"), width=2, align="center"),
            dbc.Col(dbc.Button("Update Options", className="sm",
                    id="update_mapping_candidates_experiments_dropdown"), width=2, align="center"),
        ], justify="start"),
        html.P(),
        html.Div(id='div_mapping_candidates_experiments_dropdown', style = {"width": "95%"}),
        
        html.Div(id='div_mapping_candidates_traits_dropdown', style = {"width": "95%"}),

        dbc.Row([
            dbc.Col([
                dbc.Label("Select which scores to include (not implemented yet):"),
                dbc.Checklist(options=[{"label": str(n), "value": n} for n in range(6)],
                    value=[5], id="candidate_score_checklist", 
                    inline=True, style = {"width": "95%"}
                ),]),
            dbc.Col([    
                dbc.Label("Options for selected candidates"),
                dbc.Checklist(
                    options=[
                        {"label": "Show manhattan plots?", "value": 1},
                        {"label": "Show clustered traits?", "value": 2},
                    ], id="switches-inline-input",
                    value=[1,2], inline=True, switch=True,
                ),]),
            dbc.Col([
                dbc.Label("Select grouped bins"),
                html.Br(),
                dbc.Button("Light", color="info", className="me-1", id="select_table_bins"),])
        ]),
        html.Div(id='div_mapping_candidates_selected_table', style = {"width": "95%"}),
    ])

    @dashapp.callback(
        Output('div_mapping_candidates_experiments_dropdown', 'children'),
        Input('update_mapping_candidates_experiments_dropdown', 'n_clicks'),)
    def update_mapping_candidates_experiments_dropdown(n_clicks):
        return(dcc.Dropdown(
                id='mapping_candidates_experiments_dropdown',
                options=distinct_db_vals(sqnce_path, "mapping_candidates", "experiment",[["All","all"]]),
                    #[["All","all"], ["Complete","complete"],["Skipped","skipped"],["Unannotated","unannotated"]]),
                value="all",
                multi=True,
                searchable=True,
                style = {"width": "95%"},
            ),)

    @dashapp.callback(
        Output('div_mapping_candidates_traits_dropdown', 'children'),
        Input('mapping_candidates_experiments_dropdown', 'value'),
        prevent_initial_call=True)
    def coordinates_genotypes_dropdown_div(selected_experiments):
        if "all" in selected_experiments:
            options = distinct_db_vals(sqnce_path, 'mapping_candidates',
                "trait", [["All", "all"]])
        elif len(selected_experiments)==0:
                return html.P("Please select an experiment.")
        else:
            options = []
            for experiment in selected_experiments:
                # Looks like it is removing duplicates in the dcc.Dropdown automatically
                options += return_column_if(sqnce_path, 'mapping_candidates', 'trait',
                    "experiment", experiment, [["All", "all"]])
        return html.Div([dcc.Dropdown(
            id='mapping_candidates_traits_dropdown',
            options=options,
            value="all",
            multi=True,
            searchable=True,
            style = {"width": "95%"},)
    ])

    @dashapp.callback(
        Output('div_mapping_candidates_selected_table', 'children'),
        Input('mapping_candidates_traits_dropdown', 'value'),
        Input('mapping_candidates_experiments_dropdown', 'value'),
        prevent_initial_call=True)
    def mapping_parser_candidates_table(selected_traits, selected_experiments):
        con = sqlite3.connect(sqnce_path)
        if "all" in selected_experiments and "all" in selected_traits:
            df = pd.read_sql_query('''SELECT * 
                                      FROM mapping_candidates''', con)  
        elif "all" not in selected_experiments and len(selected_experiments)>0:
            df = pd.DataFrame()
            for experiment in selected_experiments:
                df = pd.concat([df, pd.read_sql_query('''SELECT * 
                                      FROM mapping_candidates
                                      WHERE experiment="{0}"'''.format(experiment), con)])
        elif len(selected_experiments) == 0:
            return(html.P("Please select an experiment to show candidate genes."))
        df["pval"] = df["pval"].apply(lambda x: np.round(-np.log10(x), 2))
        df["bin"] = df["pos"].apply(lambda x: round(x, -5)/100000)
        df = df.drop_duplicates() # Just in case
        df = df.sort_values(["chrom", "pos"])
        # https://dash.plotly.com/datatable/tooltips
        tooltip_data=[{
            column: {'value': '', 'type': 'markdown'}
            for column, value in row.items()
        } for row in df.to_dict('records')],
        count = 0
        for row in df["trait"]:
            cursorObj = con.cursor()
            cursorObj.execute('''SELECT feature FROM "mapping_clusters" WHERE trait=?''', tuple([row]))
            selected = cursorObj.fetchall()
            if selected == []:
                count += 1
                continue
            else:
                tooltip_data[0][count]["trait"]["value"] = " ".join([x[0] for x in selected])
            count += 1
            
        # Tooltips are not correctly placed if DataTable does the sorting unfortunately

        # Not sure why it saves it as a tuple but that's why I'm using the [0]
        columns = [{"name": i, "id": i, "deletable": False, "selectable": False, "presentation": "input"} for i in df.columns],
        #print(columns)
        columns[0][-1]["presentation"] = "markdown"
        return dash_table.DataTable(
                    id='selected_candidate_table',
                    columns=columns[0],
                    # https://dash.plotly.com/datatable/reference
                    #sort_by = [{"column_id": "chrom", "direction": "asc"},
                    #            {"column_id": "pos", "direction": "asc"}],
                    data=df.to_dict('records'),
                    tooltip_data=tooltip_data[0],
                    css=[{
                        'selector': '.dash-table-tooltip',
                        'rule': 'background-color: black; font-family: monospace; color: white'}, 
                        {"selector": ".show-hide", "rule": "display: none"}], # removes "Toggle Columns" button 
                    markdown_options= {"html": True},
                    editable=False,
                    #filter_action="native",
                    sort_action="native",
                    sort_mode="multi",
                    #column_selectable="single",
                    row_selectable="multi",
                    row_deletable=False,
                    #selected_columns=[],
                    selected_rows=[],
                    page_action="native",
                    page_current= 0,
                    page_size= len(df),#100,
                    hidden_columns=["effect", "gene_end", "snp"], #"effect","query", "trait", "experiment","chrom","pos"], #"experiment", "trait",
                )

    @dashapp.callback(
        Output('selected_candidate_table', 'selected_rows'),
        Input('select_table_bins', 'n_clicks'),
        [State('selected_candidate_table', 'selected_rows'),
        State('selected_candidate_table', 'data')], 
        prevent_initial_call=True)
    def mapping_parser_candidates_table(n_clicks, selected_rows, rows):
        df = pd.DataFrame.from_dict(rows)
        chroms = df.iloc[selected_rows, :]["chrom"].values.tolist()
        bins = df.iloc[selected_rows, :]["bin"].values.tolist()
        selection = set()
        #print(chroms, bins)
        for i in range(len(bins)):
            #print("PRINT", chroms[i], bins[i])
            matching_ix = df[(df["chrom"]==chroms[i]) & (df["bin"]==bins[i])].index
            #print(matching_ix)
            selection = selection.union(set(matching_ix))
        selection = list(selection)
        return(selection)
        
    ###########################################################################
    ###########################################################################

    tab_mapping_summary_info = html.Div([
        dbc.Label("This tab will help access all the mapping candidates."),
        # To start I can do a dropboxwith "all" and then groupby and count
        # And otherwise sort by distinct symbols, gene IDs
        dbc.Row([
            dbc.Col(html.P("Select which experiments to show:"), width=2, align="center"),
            dbc.Col(dbc.Button("Update Options", className="sm",
                    id="update_mapping_candidates_experiments_dropdown"), width=2, align="center"),
        ], justify="start"),
        html.P(),
        html.Div(id='div_mapping_candidates_experiments_dropdown', style = {"width": "95%"}),
        
        html.Div(id='div_mapping_candidates_single_traits_dropdown', style = {"width": "95%"}),

        dbc.Row([
            dbc.Col([
                dbc.Label("Select which scores to include (not implemented yet):"),
                dbc.Checklist(options=[{"label": str(n), "value": n} for n in range(6)],
                    value=[5], id="candidate_score_checklist", 
                    inline=True, style = {"width": "95%"}
                ),]),
            dbc.Col([    
                dbc.Label("Options for selected candidates"),
                dbc.Checklist(
                    options=[
                        {"label": "Show manhattan plots?", "value": 1},
                        {"label": "Show clustered traits?", "value": 2},
                    ], id="switches-inline-input",
                    value=[1,2], inline=True, switch=True,
                ),]),
            dbc.Col([
                dbc.Label("Select grouped bins"),
                html.Br(),
                dbc.Button("Light", color="info", className="me-1", id="select_table_bins"),])
        ]),
        html.Div(id='div_mapping_candidates_selected_table', style = {"width": "95%"}),
    ])

    @dashapp.callback(
        Output('div_mapping_candidates_single_traits_dropdown', 'children'),
        Input('mapping_candidates_experiments_dropdown', 'value'),
        prevent_initial_call=True)
    def coordinates_genotypes_dropdown_div(selected_experiments):
        if "all" in selected_experiments:
            options = distinct_db_vals(sqnce_path, 'mapping_candidates',
                "trait")
        elif len(selected_experiments)==0:
                return html.P("Please select an experiment.")
        else:
            options = []
            for experiment in selected_experiments:
                # Looks like it is removing duplicates in the dcc.Dropdown automatically
                options += return_column_if(sqnce_path, 'mapping_candidates', 'trait',
                    "experiment", experiment)
        return html.Div([dcc.Dropdown(
            id='mapping_candidates_traits_dropdown',
            options=options,
            #value="all",
            multi=False,
            searchable=True,
            style = {"width": "95%"},)
    ])

    ###########################################################################
    ###########################################################################



    # Query to find neighboring genes
    def get_SNP_neighbors(genotype, chromosome, coordinate, distance):
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
                        '''.format(genotype, chromosome, coordinate-distance, coordinate+distance), con)
        # Should check why it returns the same row twice, probably need to correct the query
        if df.empty:
            return df
        df = df.drop_duplicates()
        # This might confusing giver the original SNP name doesn't necessarily match
        # df.insert(0, 'Query', pd.Series(["_".join([chromsome, str(coordinate)]) for x in range(len(df.index))]))
        df = df.drop("gene_chr", axis=1)
        df.insert(0, 'chrom', pd.Series([chromosome for x in range(len(df.index))]))
        df.insert(0, 'pos', pd.Series([coordinate for x in range(len(df.index))]))
        def calculate_distance(row):
            left = min([row["gene_start"], row["gene_end"]])
            right = max([row["gene_start"], row["gene_end"]])
            if left <= row["pos"] <= right:
                return 0
            else:
                return min([abs(left-row["pos"]), abs(row["pos"]-right)])
        df["distance"] = df.apply(calculate_distance, axis=1)
        return df

    def multiple_row_select(con, table, col_name, entity_list):
        ls = []
        for entity in entity_list:
            cursorObj = con.cursor()
            cursorObj.execute('''SELECT gene_id, {1} 
                                FROM {0}
                                WHERE gene_id =  ?  '''.format(table, col_name), (entity,))
            # (name,) - need the comma to treat it as a single item and not list of letters
            selected = cursorObj.fetchall()
            if selected == []:
                ls.append("")
            else:
                ls.append(selected[0][1])    
        return(ls)

   
