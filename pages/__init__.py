import dash
from flask import Flask, render_template, redirect, url_for, request
from flask.helpers import get_root_path

from flask import Flask, send_from_directory, session

from flask_login import LoginManager, login_user
from flask_wtf import FlaskForm
from wtforms import StringField, PasswordField, SubmitField

from urllib.parse import urlparse, urljoin

import dash_bootstrap_components as dbc
import dash
import pandas as pd 

from config import BaseConfig

import os

class LoginForm(FlaskForm):
    username = StringField('Username')
    password = PasswordField('Password')
    submit = SubmitField('Submit')

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

    #login_manager = LoginManager()
    #@login_manager.user_loader
    #def load_user(user_id):
    #    return User.get(user_id)

    #def is_safe_url(target):
    #    ref_url = urlparse(request.host_url)
    #    test_url = urlparse(urljoin(request.host_url, target))
    #    return test_url.scheme in ('http', 'https') and \
    #        ref_url.netloc == test_url.netloc


    # Route for handling the login page logic
    @server.route('/login', methods=['GET', 'POST'])
    def login():
        error = None
        if request.method == 'POST':
            if request.form['username'] != 'admin' or request.form['password'] != 'admin':
                error = 'Invalid Credentials. Please try again.'
            else:
                session['username'] = "admin"
                return redirect(url_for('/'))
        return render_template('login.html', error=error)

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