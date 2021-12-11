import os
basedir = os.path.abspath(os.path.dirname(__file__))
print("Base dir is:", basedir)

class BaseConfig:
    SQLALCHEMY_DATABASE_URI = os.environ.get('DATABASE_URL')
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    SECRET_KEY = 'SECRET_KEY'
