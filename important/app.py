"pip install --upgrade dash dash-core-components dash-html-components dash-renderer dash_bootstrap_components"

"""
A simple app demonstrating how to dynamically render tab content containing
dcc.Graph components to ensure graphs get sized correctly. We also show how
dcc.Store can be used to cache the results of an expensive graph generation
process so that switching tabs is fast.
"""
#con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db')
# con = sqlite3.connect('SQNce.db')

import time

import pandas as pd
from dash import dash_table as dt
from dash import dash_table

from flask import Flask
import sqlite3
import dash
import dash_bootstrap_components as dbc
from dash import dcc
from dash import html
import numpy as np
import plotly.graph_objs as go
from dash.dependencies import Input, Output

# Normally, Dash creates its own Flask server internally. By creating our own,
# we can create a route for downloading files directly:



server = Flask(__name__)
app = dash.Dash(server=server, external_stylesheets=[dbc.themes.BOOTSTRAP])

# Load callbacks
# https://community.plotly.com/t/splitting-callback-definitions-in-multiple-files/10583/2
# import call_layout
# from call_layout import *

#from call_annotations import *

if __name__ == "__main__":
    app.run_server(debug=True, host='127.0.0.1', port=8881)
#    app.run_server()