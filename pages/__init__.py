import dash
from flask import Flask
from flask.helpers import get_root_path

from flask import Flask, send_from_directory
import dash_bootstrap_components as dbc
import dash
import pandas as pd 

from dash_router import Router
from config import BaseConfig

import os

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

    from pages.Apps.layout import layout as apps_layout
    from pages.Apps.callbacks import register_callbacks as register_app_callbacks
    register_dashapp(server, 'apps', 'apps', apps_layout, register_app_callbacks)

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

    # Can't get this router function to work
    if title == "apps":
        print("Trying to register router")
        router = Router()
        router.register_callbacks(my_dashapp)
        
        if "plantapp" not in os.getcwd().lower():
            UPLOAD_DIRECTORY = "/home/eporetsky/plantapp/pages/Apps/download" # for server hosting
        else:
            UPLOAD_DIRECTORY = os.path.join(os.getcwd(), "pages", "Apps", "download") # for personal computer
        print("Router upload dir:", UPLOAD_DIRECTORY)
        @router.route('/apps/download/<path:path>')
        def download(path):
            """Serve a file from the upload directory."""
            print("If triggered path is:", path)
            return send_from_directory(UPLOAD_DIRECTORY, path, as_attachment=True)


    with app.app_context():
        my_dashapp.title = title
        my_dashapp.layout = layout
        register_callbacks_fun(my_dashapp)