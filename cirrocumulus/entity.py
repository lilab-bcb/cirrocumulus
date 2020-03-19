class Entity:
    def __init__(self, id, d):
        self.id = id
        self.d = d

    def __contains__(self, item):
        return item in self.d

    def __getitem__(self, item):
        return self.d[item]

    def get(self, key, default):
        return self.d.get(key, default)
