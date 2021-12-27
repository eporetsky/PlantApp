import sqlite3
import os 
import dash_bootstrap_components as dbc
from dash import Dash, dcc, html, Input, Output, State, dash_table




def register_callbacks(dashapp):
    from dash.dependencies import Input, Output

    import pandas as pd
    from collections import OrderedDict
    import dash_bootstrap_components as dbc
    from dash import Dash, dcc, html, Input, Output, State, dash_table

    import csv

    from flask import Flask
    import sqlite3
    import numpy as np
    import plotly.graph_objs as go

    import os
    import zlib