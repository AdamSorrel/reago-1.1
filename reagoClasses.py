import json

###############################################################################
# Classes definition
###################################################################################

class JSONObject:
    # Initialising empty dictionary
    def __init__(self):
        self.__dict__ = {}
    def update(self, d):
    # If it exist, and old key is removed and replaced by the new one
        for key in d:
            if key in self.__dict__:
                del self.__dict__[key]
        d.update(self.__dict__)
        self.__dict__ = d
    def read(self, input):
    # Reading a json file.
        with open(input, 'r') as fj:
            d = json.load(fj)
            d.update(self.__dict__)
            self.__dict__ = d
    def write(self, output):
    # Writing a json file
        with open(output, 'w') as fj:
            json.dump(self.__dict__, fj, ensure_ascii=False,separators=(',\n', ':'))

class databases:
    def __init__(self):
    # Initialising a list of database names
        self.__dict__ = {}
    #def dt(self):
        # Setting up a special dtype
        #self.dt = h5py.special_dtype(vlen=bytes)
    def newGroup(self, name, db):
        # If it exist, and old key is removed and replaced by the new one
        if name in self.__dict__:
                del self.__dict__[name]
    def add(self, name, seq):
        length = len(seq)
        d = {name : db.create_dataset(name, shape=(100,), dtype=h5py.special_dtype(vlen=bytes))}
        d.update(self.__dict__)
        self.__dict__ = d
