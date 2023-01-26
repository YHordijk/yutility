class Config:
    def __init__(self, path):
        self.data = {}
        self.path = path
        self.read()

    def read(self):
        with open(self.path) as conf:
            lines = [line.strip() for line in conf.readlines()]

        for line in lines:
            self.data[line.split('=')[0]] = line.split('=')[1]

    def save(self):
        with open(self.path, 'w+') as conf:
            for key, val in self.data.items():
                conf.write(f'{key} = {val}\n')


    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        self.data[key] = value
