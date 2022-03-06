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
    import uuid

    from collections import OrderedDict

    import dash
    from dash import no_update # https://community.plotly.com/t/error-expected-the-output-type-to-be-a-list-or-tuple-but-got-none/34744/6
    from flask import Flask, send_from_directory
    from urllib.parse import quote as urlquote

    from io import StringIO
    from Bio import SeqIO
    from Bio import Phylo
    from Bio import AlignIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    # import numpy as np
    from io import BytesIO
    import base64

    from flask_login import current_user
    from flask import session

    # import pandas as pd
    from plotly.tools import mpl_to_plotly
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
        print(session.get('username', None))
        pathname = pathname.split("/")[-1]
        if pathname == "simple_tree":
            return tab_simple_tree
        if pathname == "genome_graph":
            return tab_genome_graph



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
        # type: (plt.Figure) -> str
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