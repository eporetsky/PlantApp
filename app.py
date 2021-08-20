
import subprocess
import os
#import gzip
import zlib

from collections import OrderedDict

import dash
from dash import no_update # https://community.plotly.com/t/error-expected-the-output-type-to-be-a-list-or-tuple-but-got-none/34744/6
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Output, State, Input

#import flask
#from flask.helpers import send_file
from flask import Flask, send_from_directory

from urllib.parse import quote as urlquote

import sqlite3
#from sqlite3 import Error

from io import StringIO
from Bio import SeqIO
from Bio import Phylo
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# import numpy as np
from io import BytesIO
import base64


# import pandas as pd
from plotly.tools import mpl_to_plotly
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


UPLOAD_DIRECTORY = "download"

#if not os.path.exists(UPLOAD_DIRECTORY):
#    os.makedirs(UPLOAD_DIRECTORY)

# Normally, Dash creates its own Flask server internally. By creating our own,
# we can create a route for downloading files directly:
server = Flask(__name__)
app = dash.Dash(server=server)


@server.route("/download/<path:path>")
def download(path):
    """Serve a file from the upload directory."""
    return send_from_directory(UPLOAD_DIRECTORY, path, as_attachment=True)


colors = {
    'background': '#111111',
    'text': '#7FDBFF'
}


app.layout = html.Div([
    
    dcc.Tabs(
        id="tabs-with-classes",
        value='tab-2',
        parent_className='custom-tabs',
        className='custom-tabs-container',
        children=[
            dcc.Tab(
                label='Tab one',
                value='tab-1',
                className='custom-tab',
                selected_className='custom-tab--selected'
            ),
            dcc.Tab(
                label='Tab two',
                value='tab-2',
                className='custom-tab',
                selected_className='custom-tab--selected'
            ),
            dcc.Tab(
                label='Tab three, multiline',
                value='tab-3', className='custom-tab',
                selected_className='custom-tab--selected'
            ),
            dcc.Tab(
                label='Tab four',
                value='tab-4',
                className='custom-tab',
                selected_className='custom-tab--selected'
            ),
        ]),
    html.Div(id='tabs-content-classes'),
    
    dcc.Input(id='gene_list', value='Paste gene list', type='text'),
    html.Button('Get gene list', id='fasta-button', type='submit'),#, children='Submit'),
    html.H2("Download Fasta File"),
    html.Ul(id="fasta-download"),
    html.Button('Prepare alignment and tree', id='aln-button', type='submit'),
    html.H2("Download Alignment and Tree"),
    html.Ul(id="aln-download"),
    html.Ul(id="tree-download"),
    html.Ul(id="aln-tree-png-download"),
    html.Div([html.Img(id = 'aln-tree-figure', src = '')], id='plot_div'),
    ])

@app.callback(Output('tabs-content-classes', 'children'),
              Input('tabs-with-classes', 'value'))
def render_content(tab):
    if tab == 'tab-1':
        return html.Div([
            html.H3('Tab content 1')
        ])
    elif tab == 'tab-2':
        return html.Div([
            html.H3('Tab content 2')
        ])
    elif tab == 'tab-3':
        return html.Div([
            html.H3('Tab content 3')
        ])
    elif tab == 'tab-4':
        return html.Div([
            html.H3('Tab content 4')
        ])

@app.callback(
    Output("fasta-download", "children"),
    [Input('fasta-button', 'n_clicks')],
    [State('gene_list', 'value')])
def fasta_update(clicks, input_value):
    if clicks is not None:
        try:
            #os.makedirs("/home/eporetsky/plantapp/worked_2")
            con = sqlite3.connect('/home/eporetsky/plantapp/PlantApp.db')
            input_value = ("Zm00019ab401890", "Zm00028ab426680", "Zm00020ab421830", "Zm00030ab419340",) #(input_value.split(" ")
            #print(clicks, input_value)
            protein_seq_select(con, input_value)
            #print("ran protein")
            con.close()
            print("closed connection")
            # uses a separate function to generate the download link
            return [html.Li(file_download_link("selected.fasta"))]
        except:
            input_value = "didn't work"

        return(input_value)

def protein_seq_select(con, entity_list):
    #os.makedirs("/home/eporetsky/plantapp/"+entity_list[0])
    od = OrderedDict()
    for entity in entity_list:
        cursorObj = con.cursor()
        cursorObj.execute('''SELECT protein_variant, protein_sequence
                             FROM protein_seqs
                             WHERE protein_variant =  ?  ''', (entity,))
        # (name,) - need the comma to treat it as a single item and not list of letters
        selected = cursorObj.fetchall()[0]
        record = SeqRecord(Seq(zlib.decompress(selected[1]).decode(encoding='UTF-8')), id=selected[0], name="", description="")
        od[selected[0]] = record

        with open("/home/eporetsky/plantapp/download/selected.fasta", 'w') as handle:
            SeqIO.write(od.values(), handle, 'fasta')

# https://stackoverflow.com/questions/54443531/downloading-dynamically-generated-files-from-a-dash-flask-app
def file_download_link(filename):
    """Create a Plotly Dash 'A' element that downloads a file from the app."""
    location = "/download/{}".format(urlquote(filename))
    #location = "/home/eporetsky/plantapp/download/selected.fasta"
    return html.A(filename, href=location)


@app.callback(
    Output("aln-download", "children"),
    Output("tree-download", "children"),
    Output("aln-tree-png-download", "children"),
    Output("aln-tree-figure", "src"),
    [Input('aln-button', 'n_clicks')])
def alignment_update(clicks):
    ctx = dash.callback_context
    if (not ctx.triggered and not ctx.triggered[0]['value'] == 0):
        return (no_update, no_update, no_update, no_update)
    if clicks is not None:
        try:
            print("running famsa")
            subprocess.run(["/home/eporetsky/plantapp/./famsa /home/eporetsky/plantapp/download/selected.fasta /home/eporetsky/plantapp/download/selected.aln", "arguments"], shell=True)
            print("running clipkit")
            subprocess.run(["python3 /home/eporetsky/plantapp/clipkit /home/eporetsky/plantapp/download/selected.aln", "arguments"], shell=True)
            print("running fasttree")
            subprocess.run(["/home/eporetsky/plantapp/./fasttree /home/eporetsky/plantapp/download/selected.aln.clipkit > /home/eporetsky/plantapp/download/selected.tree", "arguments"], shell=True)
            ########################################
            # Now we should have a tree and an alignment file in the download folder
            # In 3.x, there is no longer a .next method; use the built-in function:
            tree = next(Phylo.parse('/home/eporetsky/plantapp/download/selected.tree', 'newick'))
            align = AlignIO.read("/home/eporetsky/plantapp/download/selected.aln", "fasta")
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(len(align), 8), constrained_layout=True)
            #fig.suptitle('Horizontally stacked subplots')
            # https://stackoverflow.com/questions/29419973/in-python-how-can-i-change-the-font-size-of-leaf-nodes-when-generating-phylogen
            ########################################
            line_width = 0.1 # sometimes if the line is not thick enought it might not appear in the %inline
            box_width = 0.5
            #add rectangle to plot
            def names(terminal):
                return terminal.name
            phylo_ordered_list = list(map(names, tree.get_terminals()))

            for aln in align:
                aln_count = len(align) - phylo_ordered_list.index(aln.id) - 1
                ax2.add_patch(Rectangle((0, aln_count - line_width/2), len(align[0].seq), line_width))
                rect_ls = alignment_to_range_list(aln.seq)
                for rect in range(0, len(rect_ls), 2):
                    ax2.add_patch(Rectangle((rect_ls[rect], aln_count-box_width/2), rect_ls[rect+1] - rect_ls[rect], box_width))
                aln_count -= 1
            ax2.set(yticklabels=[])
            ax2.set(ylabel=None)
            ax2.tick_params(left=False)  # remove the ticks
            ax2.plot(0)
            ########################################
            # Keep this reverse order because I Phylo.draw uses plt.show() (I assume) which prevents ax2 from showing up
            matplotlib.rc('font', size=8)
            #plt.figure(figsize=(20, 6))
            Phylo.draw(tree, axes=ax1, do_show=False)
            #ax.add_patch(Rectangle((box_width, box_width), 2, 6))
            fig.savefig("/home/eporetsky/plantapp/download/selected_tree_aln.png")
            ########################################
            aln_file = [html.Li(file_download_link("selected.aln"))]
            tree_file = [html.Li(file_download_link("selected.tree"))]
            figure_file = [html.Li(file_download_link("selected_tree_aln.png"))]
            return aln_file, tree_file, figure_file, fig_to_uri(fig)
        except:
            input_value = "Make sure a fasta file exists"

        return(input_value)


def alignment_to_range_list(alignment):
    # function takes a single sequence from an alignment
    # returns a 1-d list of paired star and end coordinates for each aligned segment
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

if __name__ == '__main__':
    app.run_server(debug=True)
    #app.run_server()
    
    