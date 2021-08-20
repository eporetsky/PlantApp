# This file contains the WSGI configuration required to serve up your
# web application at http://<your-username>.pythonanywhere.com/
# It works by setting the variable 'application' to a WSGI handler of some
# description.
#
# The below has been auto-generated for your Flask project

import sys

# add your project directory to the sys.path
project_home = '/home/eporetsky/plantapp'
if project_home not in sys.path:
    sys.path = [project_home] + sys.path

# import flask app but need to call it "application" for WSGI to work
# from app import app as application  # noqa

# https://help.pythonanywhere.com/pages/DashWSGIConfig
# Dash apps have a app.server that you need to pass to WSGI
from app import server as application

#from app import app
#application = app.server


# application = app.server
