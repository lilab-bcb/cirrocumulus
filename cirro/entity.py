class Entity:
    def __init__(self, id, d):
        self.id = id
        self.d = d

    def __getitem__(self, item):
        return self.d[item]
