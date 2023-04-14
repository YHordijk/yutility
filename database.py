import sqlite3 as sql
from yutility import log, ensure_list, dictfunc
import os
import numpy as np


python_to_sql_types = {
    str: 'TEXT',
    int: 'INTEGER',
    float: 'REAL',
    bool: 'BOOL',
    None: 'NULL',
    type: 'BLOB'
}

sql_to_python_types = {
    'TEXT': str,
    'INTEGER': int,
    'REAL': float,
    'NULL': None,
    'BOOL': bool
}


class DBSelectResult:
    def __init__(self, data, columns):
        self.data = data
        self.columns = columns
        self.__iter_counter = 0

    def __getitem__(self, key):
        if isinstance(key, int):
            return self.data[key]
        
        col_idx = self.columns.index(key)
        return np.array([datum[col_idx] for datum in self.data])

    def __setitem__(self, key, values):
        if isinstance(key, int):
            # assert len(self.data[key]) == len(values), f'Cannot set values with length {len(self.data[key])} at index {key} with length {len(values)}'
            self.data[key] = values
        
        col_idx = self.columns.index(key)
        # assert len(self.data[key]) == len(values), f'Cannot set values with length {len(self.data[key])} at column {key} with length {len(values)}'
        newdata = []
        for val, row in zip(values, self.data):
            row = list(row)
            row[col_idx] = val
            newdata.append(row)
        self.data = newdata
            
    def groupby(self, mask_key: str):
        mask_data = self[mask_key]
        unqs = set(mask_data)
        for unq in unqs:
            yield DBSelectResult([x for i, x in enumerate(self) if mask_data[i] == unq], self.columns)

    def __iter__(self):
        return iter(self.data)

    def __add__(self, other):
        assert isinstance(other, DBSelectResult), f'Added object must be of type DBSelectResult, not {type(other)}'
        self.data.extend(other.data)
        return self

    def __radd__(self, other):
        assert isinstance(other, DBSelectResult), f'Added object must be of type DBSelectResult, not {type(other)}'
        self.data.extend(other.data)
        return self

    def __len__(self):
        return len(self.data)
        

class DataBase:
    def __init__(self, db_path):
        self.db_path = db_path
        self.name = os.path.split(db_path)[1]
        self.connection = sql.connect(self.db_path)
        self.cursor = self.connection.cursor()

    def __str__(self):
        tables = [f'{name}({self.get_table_size(name)}x{len(self.get_column_names(name))})' for name in self.get_table_names()]
        return f'DataBase(file={self.name}, tables=[{", ".join(tables)}])'

    def __enter__(self):
        return self

    def close(self):
        self.connection.commit()
        self.connection.close()

    def __exit__(self, *args):
        self.close()

    def __del__(self):
        try:
            self.close()
        except:
            pass

    def execute(self, command):
        while True:
            try:
                self.cursor.execute(command)
                return
            except sql.OperationalError as e:
                msg = e.args[0]
                if msg == 'database is locked':
                    log.log('Warning: database locked, retrying ...')
                else:
                    print(command)
                    raise

    def fetchall(self):
        return self.cursor.fetchall()

    def fetchone(self):
        return self.cursor.fetchone()

    def parse_type(self, typ):
        if isinstance(typ, str):
            return typ
        if isinstance(typ, np.generic):
            typ = typ.item()
        else:
            return python_to_sql_types[typ]

    def make_table(self, table_name, columns=None, types=None, primary_key=None):
        columns = columns or ['id']
        types = types or [str]

        tolist = []
        for colname, coltype in zip(columns, types):
            s = f'{colname} {self.parse_type(coltype)}'
            if colname == primary_key:
                s += ' PRIMARY KEY'
            tolist.append(s)
        cols = ', '.join(tolist)
        command = f"CREATE TABLE IF NOT EXISTS {table_name} ({cols})"
        self.execute(command)

    def drop_table(self, table_name):
        command = f'DROP TABLE IF EXISTS {table_name}'
        self.execute(command)

    def insert(self, table_name, values):
        vals = ', '.join([repr(x) for x in values])
        command = f'INSERT INTO {table_name}\nVALUES ({vals})'
        self.execute(command)

    def insert_dict(self, table_name, d, ensure_columns=False):
        d = dictfunc.remove_empty(d)
        column_names = d.keys()
        values = d.values()

        for column, value in zip(column_names, values):
            if ensure_columns:
                if isinstance(value, np.generic):
                    typ = type(value.item())
                else:
                    typ = type(value)
                self.ensure_column(table_name, column, typ)
            if column not in self.get_column_names(table_name):
                raise KeyError(f"Column '{column}' not present in table '{table_name}'. Consider calling with 'ensure_columns=True'.")

        cols = ', '.join([repr(x) for x in column_names])
        vals = ', '.join([repr(x) for x in values])

        if not cols or not vals:
            return

        command = f'INSERT INTO {table_name} ({cols})\nVALUES ({vals})'
        self.execute(command)

    def delete_row(self, table_name, where):
        self.execute(f'DELETE FROM {table_name} WHERE {where}')

    def get_table_size(self, table_name):
        self.execute(f'SELECT Count(*) FROM {table_name}')
        return self.fetchone()[0]

    def get_column_names(self, table_name):
        self.execute(f'PRAGMA table_info({table_name})')
        return [row[1] for row in self.fetchall()]

    def get_column_types(self, table_name):
        self.execute(f'PRAGMA table_info({table_name})')
        return [sql_to_python_types[row[2].upper()] for row in self.fetchall()]

    def get_table_names(self):
        self.execute("SELECT name FROM sqlite_master WHERE type='table';")
        return [name[0] for name in self.fetchall()]

    def add_column(self, table_name, column, typ):
        command = f'ALTER TABLE {table_name} ADD COLUMN {column} {self.parse_type(typ)}'
        self.execute(command)

    def delete_column(self, table_name, column):
        ...

    def ensure_column(self, table_name, column, typ):
        if column not in self.get_column_names(table_name):
            self.add_column(table_name, column, typ)

    def select(self, table_name, columns=None, where=None):
        cols = '*'
        if columns is not None:
            columns = ensure_list(columns)
            cols = ', '.join(columns)
        if cols == '*':
            columns = self.get_column_names(table_name)
        if where is not None:
            where = 'WHERE ' + where
        command = f'SELECT {cols}\n\tFROM {table_name}\n\t{where}'
        self.execute(command)
        result = self.fetchall()
        return DBSelectResult(result, columns)

    def delete_duplicates(self, table_name, columns=None):
        cols = '*'
        if columns is not None:
            columns = ensure_list(columns)
            cols = ', '.join(columns)
        command = rf"""
        DELETE FROM {table_name}
        WHERE rowid NOT IN (
          SELECT 
            MIN(rowid) 
          FROM 
            {table_name} 
          GROUP BY
            {', '.join(cols)}
        )"""
        print(command)


def merge_databases(databases, new_name):
    with DataBase(new_name) as db:
        for db_old in databases:
            for table in db_old.get_table_names():
                db.make_table(table)
                cols = db_old.get_column_names(table)
                data = db_old.select(table, '*')
                for datum in data:
                    datum_dict = {col: value for col, value in zip(cols, datum)}
                    db.insert_dict(table, datum_dict, ensure_columns=True)

    return DataBase(new_name)

