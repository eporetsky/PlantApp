import dash_core_components as dcc
import dash_html_components as html
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

search_bar = [ 
    dbc.Row([
    dbc.Col(dbc.Input(type="search", placeholder="Search")),
        dbc.Col(
            dbc.Button("Search", color="primary", className="ml-2", n_clicks=0),
            width="auto",),
        ],
    no_gutters=True,
    className="ml-auto flex-nowrap mt-3 mt-md-0",
    align="center",
    ), 
]

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
                align="center",# no_gutters=True,
            ),
            href="https://www.plantapp.org",
        ),

        dbc.NavbarToggler(id="navbar-toggler", n_clicks=0),
        dbc.Collapse(search_bar, id="navbar-collapse", navbar=True, is_open=False),
    ],
    color="dark", dark=True,
)

#table = dbc.Table.from_dataframe(df, striped=True, bordered=True, hover=True)

layout = dbc.Container(
    [
        dcc.Store(id="store"),

        #html.H1("Dynamically rendered tab content"),
        #html.Hr(),

        navbar,

        dbc.Tabs(
            [
                dbc.Tab(label="Available DBs", tab_id="available_dbs"),
                dbc.Tab(label="Annotations", tab_id="annotations"),
                dbc.Tab(label="Sequences", tab_id="proteins"),
                dbc.Tab(label="Promoters", tab_id="promoters"),
                dbc.Tab(label="Omics", tab_id="omics"),
            ],
            id="tabs",
            active_tab="scatter",
        ),
        
        html.Div(id="tab_content"),

        
      #dbc.Row([table])
    ],
    fluid=True
)