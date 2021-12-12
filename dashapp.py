from app import create_app
import os, sqlite3

server = create_app()

if __name__ == "__main__":
    server.run(debug=True, host='127.0.0.1', port=8881)