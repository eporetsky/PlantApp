
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