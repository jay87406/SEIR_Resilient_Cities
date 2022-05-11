from flask import Flask
app = Flask(__name__)

@app.route('/')
def hello_world():
    return 'hello world'

@app.route('/login')
def login():
    return 'Login'


if __name__ == '__main__':
    app.run()