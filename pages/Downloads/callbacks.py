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
    from scipy.stats import  hypergeom
    from statsmodels.stats.multitest import multipletests

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
        if pathname == "plantapp":
            return tab_downloads_plantapp

    tab_downloads_plantapp = html.Div([
        html.P("Select which data type to download the raw data from:"),
        dcc.Dropdown(
            id='downloads_select_folder',
            options=[{'label': 'Protein sequences', 'value': 'proteins'},
                     {'label': 'Promoter sequences', 'value': 'promoters'},
                     {'label': 'Gene family annotations', 'value': 'families'},
                     {'label': 'Best-Blast-Hits', 'value': 'BBHs'},
                     {'label': 'Gene Symbols', 'value': 'symbols'}],
            value="Select download type"
        ),
        html.Div(id='download_plantapp_file_div', style = {"width": "95%"}),
        dcc.Download(id='download_plantapp_file_download')
    ])

    @dashapp.callback(Output("download_plantapp_file_div", "children"),
                      Input("downloads_select_folder", "value"),
                      prevent_initial_call=True)
    def update_output(name):
        # Couldn't get it to work with html.Li links. Needs a dash download call back instead
        # https://docs.faculty.ai/user-guide/apps/examples/dash_file_upload_download.html
        # https://towardsdatascience.com/creating-interactive-data-tables-in-plotly-dash-6d371de0942b
        # https://dash.plotly.com/dash-core-components/download
        
        fldr = os.path.join(cwd, "downloads", name)
        return(
            dash_table.DataTable(
            id='downloads_plantapp_data_table',
            columns=[{'name': 'File Name', 'id': "fl_name"}],
            data=[{'fl_name': fl} for fl in os.listdir(fldr)]
        ),
    )

    @dashapp.callback(Output("download_plantapp_file_download", "data"),
                      Input("downloads_plantapp_data_table", "active_cell"),
                      Input("downloads_plantapp_data_table", "data"),
                      State("downloads_select_folder", "value"),
                      prevent_initial_call=True)
    def func(active_cell, data, name):
        if active_cell:
            col = active_cell['column_id']
            row = active_cell['row']
            cellData = data[row][col]
            print(os.path.join("downloads", name, cellData))
            return(dcc.send_file(os.path.join(cwd, "downloads", name, cellData)))