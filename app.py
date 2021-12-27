from pages import create_app
import os, sqlite3

server = create_app()

if __name__ == "__main__":
    server.run(debug=True, host='127.0.0.1', port=8881)

# Search and replace when deploying on server
# con = sqlite3.connect(os.getcwd()+'\SQNce.db') # for windows
# con = sqlite3.connect(os.getcwd()+'/SQNce.db') # on linux
# con = sqlite3.connect('/home/eporetsky/plantapp/SQNce.db') # deploy with this