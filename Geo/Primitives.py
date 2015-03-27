'''
Created on Jan 29, 2015

@author: Kutach
'''

class GeoOb(object):
    '''
    This is the superclass for all geometric objects processed by the proof system.
    '''

    def __init__(self, label=None, *params, **keywordargs):
        '''
       Initializes the directory that will hold this objects properties.
        '''
        self.label = label
        self.tags = keywordargs.copy()

    def __repr__(self):
        s = self.label + " { "
        for key, value in self.tags.items():
            s += str(key) + " => " + str(value) + ", "
        s += " } "
        return s

class Spot(GeoOb):
    '''
    This is the superclass for all geometric objects that are localized within a
    point-sized infinitessimal region. It includes points, empty nested locations,
    tanginos, circlinos, infinitessimal planar regions (beams), and possibly more.
    '''
    def __init__(self, *params, **keywordargs):
        GeoOb.__init__(self, *params, **keywordargs)

if __name__ == '__main__':
    #  Test code
    X = Spot('A', hasPoint=True, hasTangino=False, hasCirculino=False)
    print(X)
    Y = Spot('B')
    print(Y)
