import dash
from flask import Flask
from flask.helpers import get_root_path

from flask import Flask
import dash_bootstrap_components as dbc
import dash
import pandas as pd 


from config import BaseConfig

def create_app():
    server = Flask(__name__)
    server.config.from_object(BaseConfig)

    from pages.Index.layout import layout as layout0
    from pages.Index.callbacks import register_callbacks as register_callbacks0
    register_dashapp(server, 'Index', '/', layout0, register_callbacks0)
    register_dashapp(server, 'Index', 'index', layout0, register_callbacks0)

    from pages.SQNce.layout import layout as layout1
    from pages.SQNce.callbacks import register_callbacks as register_callbacks1
    register_dashapp(server, 'SQNce', 'SQNce', layout1, register_callbacks1)

    #register_dashapp(server, 'Index', '/', layout1, register_callbacks1)

    
    return server


def register_dashapp(app, title, base_pathname, layout, register_callbacks_fun):
    # Meta tags for viewport responsiveness
    meta_viewport = {"name": "viewport", "content": "width=device-width, initial-scale=1, shrink-to-fit=no"}

    external_stylesheets = [dbc.themes.BOOTSTRAP]

    if base_pathname != "/":
        url_base_pathname = f'/{base_pathname}/' 
    else:
        url_base_pathname = "/"

    my_dashapp = dash.Dash(__name__,
                           server=app,
                           url_base_pathname=url_base_pathname,
                           external_stylesheets=external_stylesheets
                           #assets_folder=get_root_path(__name__) + f'/{base_pathname}/assets/',
                           #meta_tags=[meta_viewport]
                           )

    with app.app_context():
        my_dashapp.title = title
        my_dashapp.layout = layout
        register_callbacks_fun(my_dashapp)