import sqlite3

class Storer:

    def __init__(self):
        self.connection = sqlite3.connect('src/db/proteomics.db')

    def __create_tables(self):
        pass

# If running this file: check if the database connection works.
if __name__ == '__main__':
    storer = Storer()
