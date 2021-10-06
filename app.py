"pip install --upgrade dash dash-core-components dash-html-components dash-renderer dash_bootstrap_components"

"""
A simple app demonstrating how to dynamically render tab content containing
dcc.Graph components to ensure graphs get sized correctly. We also show how
dcc.Store can be used to cache the results of an expensive graph generation
process so that switching tabs is fast.
"""
import time

import pandas as pd
import dash_table as dt

from flask import Flask
import sqlite3
import dash
from dash import dash_table
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

import call_blastree
from call_variables import search_bar

PLOTLY_LOGO = "https://images.plot.ly/logo/new-branding/plotly-logomark.png"

navbar = dbc.Navbar(

    children = [
        html.A(
            # Use row and col to control vertical alignment of logo / brand
            dbc.Row(
                [dbc.Col(html.Img(src=PLOTLY_LOGO, height="30px")),
                 dbc.Col(dbc.NavbarBrand("PlantApp", className="ml-2")),
                 dbc.Col(dbc.NavLink("About", href="/about", className="ml-2")),
                 dbc.Col(dbc.NavLink("SQNce", href="/sqnce", className="ml-2")),
                 dbc.Col(dbc.NavLink("Queries", href="/queries", className="ml-2")),
                 dbc.Col(dbc.NavLink("Tools", href="/tools", className="ml-2")),
                 dbc.Col(dbc.NavLink("Downloads", href="/downloads", className="ml-2")),
                 dbc.Col(dbc.NavLink("Feedback", href="/feedback", className="ml-2")),
                 ],
                align="center", no_gutters=True,
            ),
            href="https://www.plantapp.org",
        ),

        dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
        dbc.Collapse(search_bar, id="navbar-collapse", navbar=True, is_open=False),
    ],
    color="dark", dark=True,
)

#table = dbc.Table.from_dataframe(df, striped=True, bordered=True, hover=True)


app.layout = dbc.Container(
    [
        dcc.Store(id="store"),

        #html.H1("Dynamically rendered tab content"),
        html.Hr(),

        navbar,

        dbc.Tabs(
            [
                dbc.Tab(label="Available DBs", tab_id="available_dbs"),
                dbc.Tab(label="Annotations", tab_id="annotations"),
                dbc.Tab(label="Sequences", tab_id="sequences"),
            ],
            id="tabs",
            active_tab="scatter",
        ),
        
        html.Div(id="tab_content"),

        
      #dbc.Row([table])
    ]
)

from call_layout import *
#from call_annotations import *


# Load callbacks
# https://community.plotly.com/t/splitting-callback-definitions-in-multiple-files/10583/2
#import call_layout

if __name__ == "__main__":
    app.run_server(debug=True, port=8880)
#    app.run_server()