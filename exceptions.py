class MoleculeNotFoundError(Exception):
    def __init__(self, path):
        self.path = path

    def __str__(self):
        return f'Molecule {self.path} was not found!'


class ReactionGeneratorError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message
