'''
Created on Mar 9, 2015

@author: 10034888
'''
import numpy as np

def signpermutation(t):
    sl = list(t).copy()
    assert len(sl) == len(set(sl))
    swaps = 0
    swapflag = True
    while swapflag:
        swapflag = False
        for i in range(0, len(t) - 1):
            assert sl[i] != sl[i + 1]
            if sl[i] > sl[i + 1]:
                swaps += 1
                swapflag = True
                sl[i], sl[i + 1] = sl[i + 1], sl[i]

    if swaps % 2 == 1:
        return -1
    else:
        return +1

def choose(s, n):
    """
    Returns a list of all the possible n-length subsets that can be formed from the unique elements of sequence s.
    """
    resultset = set([])
    if n > 0 and len(s) >= n:
        recursiveset = set(s)
        member = recursiveset.pop()
        if n == 1:
            resultset.add(frozenset([member]))
        else:
            for memberset in choose(recursiveset, n - 1):
                newset = set(memberset)
                newset.add(member)
                resultset.add(frozenset(newset))
        resultset.update(choose(recursiveset, n))
    return resultset

def flatten(l, ltypes=(list, tuple)):
    """
    Flattens the iterable l (in place if mutable)
    For example, flatten(["first",["second","third"],[["fourth"]]])
    returns ["first", "second", "third", "fourth"].
    Lifted from http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
    """
    ltype = type(l)
    l = list(l)
    i = 0
    while i < len(l):
        while isinstance(l[i], ltypes):
            if not l[i]:
                l.pop(i)
                i -= 1
                break
            else:
                l[i:i + 1] = l[i]
        i += 1
    return ltype(l)

def multisign(s):
    """For an input sequence of length n, this function returns the 2^(n-1) various sums corresponding to all 
    the possible signs (positive or negative) each member could have.     
    """
    result = []
    if len(s) > 0:
        if len(s) == 1:
            return [s[0]]
        else:
            firstnum = s[0]
            r = multisign(s[1:])
            result = flatten([[firstnum + num, -firstnum + num] for num in r])
    return result

class Bracketpoly(object):
    """
    Provides a general representation of a multihomogenous bracket polynomial in a way that minimizes
    illusory degrees of freedom. For example, the polynomial in the Grassman-Pluecker relation [ab][cd]-[ac][bd]+[ad][bc] = 0
    is represented as dimension 2, degree 2, and pattern "f???" where any letter (here f) stands for a single point that is held the same
    in all terms, and ? stands for one of the other variables that varies from term to term. This representation automatically
    incorporates the facts that (1) all the variables are distinct and (1) that switching two brackets does nothing and (3) that 
    switching two adjacent variables in a single bracket results in a change of sign. The default pattern, None, is equivalent to a
    string of ?'s of length dimension*degree.
    """
    _alphabet = "abcdefghijklmnopqrstuvwxyz"

    def __init__(self, dim, deg, pattern=None):
        '''
       Initializes this multihomogenous bracket polynomial to use brackets of dimension dim, monomials of degree deg,
       and 
        '''
        self.dimension = dim
        self.degree = deg
        self.pattern = tuple([["?" * dim] for i in range(0, deg)])
        self.internalpattern = {}
        self.constantpoints = set([])
        self.fixedvariables = 0
        if pattern:
            i = 0
            for c in pattern:
                if i >= len(self.pattern):
                    print("Warning: Bracketpoly initialized with an overly long string.")
                    break
                if c != "?" and c.isalpha():
                    self.pattern[i:i + 1] = c
                    self.constantpoints.add(c)
                    self.internalpattern[self.fixedvariables] = dim
                    self.fixedvariables += 1
                i += 1
        self.constantpoints = frozenset(self.constantpoints)
        self.variablepoints = dim * deg - len(self.constantpoints)
        assert self.variablepoints >= 0
        assert self.variablepoints <= dim * deg

    def __repr__(self):
        return self.pattern

    def __print__(self, labelstring=_alphabet):
        s = ''
        variableletterindex = 0
        constantletterindex = 0
        patternindex = 0
        for deg in range(0, self.degree):
            s += '['
            for dim in range(0, self.dimension):
                letter = self.pattern[patternindex]
                patternindex += 1
                if letter != '?':
                    s += letter
                else:
                    s += labelstring[variableletterindex]
                    variableletterindex += 1
                patternindex += 1
            s += ']'
            if deg < self.degree:
                if sign > 0:
                    s += '+'
                else:
                    s += '-'
        return s

    def compute(self, np_array_vectors):
        count = self.dimension * self.degree
        total = 0
        return total

if __name__ == '__main__':

    b = Bracketpoly(2, 2, "a???")
    print(b)

#      print(multisign([2, 3]))
#      print(multisign([2, 3, 5]))

#      print(choose([], 2))
#      print(choose([1], 1))
#      print(choose([1, 2], 1))
#      print(choose([1, 2], 2))
#      print(choose([1, 2, 3, 4], 1))
#      print(choose([1, 2, 3, 4], 2))
#      print(choose([1, 2, 3, 4], 3))
#      choice = choose([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 4)
#      print(choice)
#      print(len(choice))


    assert signpermutation((0, 1, 2, 3, 4, 5)) == 1
    assert signpermutation((1, 2, 3, 4, 5, 0)) == -1
    assert signpermutation((1, 3, 2, 4, 5, 0)) == 1
    assert signpermutation((1, 2, 3, 4, 0, 5)) == 1
    assert signpermutation((1, 3, 2, 4, 0, 5)) == -1

    print("Calculating brackets for Grassmann-Plucker relations.")
    #  Construct 6 random 3d vectors.
    g = np.random.random((6, 3))
    #  Calculate 10 determinant products
    all_brackets = {}
    all_terms = {}
    label = {}
    togethertotals = np.zeros((6, 6))
    separatetotals = np.zeros((6, 6))
    sixpoints = (0, 1, 2, 3, 4, 5)
    septuples = {}
    for i in sixpoints:
        for j in sixpoints:
            septuples[i * 6 + j] = tuple()

    for i in range(0, 4):
        for j in range(i + 1, 5):
            for k in range(j + 1, 6):
                #  Uniquely identify the triplet using a key made from prime factors
                keys = (i, j, k)
                skeys = sorted([i, j, k])
                key = 2 ** (skeys[0] + 1) * 3 ** (skeys[1] + 1) * 5 ** (skeys[2] + 1)
                #  Construct the key for the matching face
                mkeys = tuple([item for item in sixpoints if item not in keys])
                smkeys = tuple(sorted(mkeys))
                mkey = 2 ** (smkeys[0] + 1) * 3 ** (smkeys[1] + 1) * 5 ** (smkeys[2] + 1)
                assert set(skeys).isdisjoint(set(smkeys))
                assert len(set(skeys).union(set(smkeys))) == 6

                #  Construct matrix
                m = np.zeros((3, 3))
                for col in (0, 1, 2):
                    for row in (0, 1, 2):
                        m[row][col] = g[keys[col]][row]
                d1 = np.linalg.det(m)
                assert d1 != 0.0
                #  Construct matrix for the matching face
                m = np.zeros((3, 3))
                for col in (0, 1, 2):
                    for row in (0, 1, 2):
                        m[row][col] = g[mkeys[col]][row]
                d2 = np.linalg.det(m)
                assert d2 != 0.0
                d = d1 * d2

                assert key not in all_brackets.keys()
                all_brackets[key] = d

                if key + mkey in all_terms:
                    #  print(d, d1, d2, signpermutation(keys + mkeys), keys, mkeys)
                    assert all_terms[key + mkey] == d
                else:
                    #  print(d, d1, d2, signpermutation(keys + mkeys), keys, mkeys)
                    all_terms[key + mkey] = d
                    label[key + mkey] = str(keys) + str(mkeys)

                    d1 = d * signpermutation(keys + mkeys)
                    d2 = d * signpermutation(mkeys + keys)

                    t1 = (d1,); t2 = (d2,)

                    togethertotals[keys[0]][keys[1]] += d1
                    togethertotals[keys[1]][keys[0]] += d1
                    togethertotals[keys[0]][keys[2]] += d1
                    togethertotals[keys[2]][keys[0]] += d1
                    togethertotals[keys[1]][keys[2]] += d1
                    togethertotals[keys[2]][keys[1]] += d1

                    togethertotals[mkeys[0]][mkeys[1]] += d2
                    togethertotals[mkeys[1]][mkeys[0]] += d2
                    togethertotals[mkeys[1]][mkeys[2]] += d2
                    togethertotals[mkeys[2]][mkeys[1]] += d2
                    togethertotals[mkeys[0]][mkeys[2]] += d2
                    togethertotals[mkeys[2]][mkeys[0]] += d2

                    septuples[keys[0] * 6 + mkeys[0]] += t1
                    septuples[keys[0] * 6 + mkeys[1]] += t1
                    septuples[keys[0] * 6 + mkeys[2]] += t1
                    septuples[keys[1] * 6 + mkeys[0]] += t1
                    septuples[keys[1] * 6 + mkeys[1]] += t1
                    septuples[keys[1] * 6 + mkeys[2]] += t1
                    septuples[keys[2] * 6 + mkeys[0]] += t1
                    septuples[keys[2] * 6 + mkeys[1]] += t1
                    septuples[keys[2] * 6 + mkeys[2]] += t1

                    septuples[mkeys[0] * 6 + keys[0]] += t2
                    septuples[mkeys[0] * 6 + keys[1]] += t2
                    septuples[mkeys[0] * 6 + keys[2]] += t2
                    septuples[mkeys[1] * 6 + keys[0]] += t2
                    septuples[mkeys[1] * 6 + keys[1]] += t2
                    septuples[mkeys[1] * 6 + keys[2]] += t2
                    septuples[mkeys[2] * 6 + keys[0]] += t2
                    septuples[mkeys[2] * 6 + keys[1]] += t2
                    septuples[mkeys[2] * 6 + keys[2]] += t2

                    separatetotals[keys[0]][mkeys[0]] += d1
                    separatetotals[keys[0]][mkeys[1]] += d1
                    separatetotals[keys[0]][mkeys[2]] += d1
                    separatetotals[keys[1]][mkeys[0]] += d1
                    separatetotals[keys[1]][mkeys[1]] += d1
                    separatetotals[keys[1]][mkeys[2]] += d1
                    separatetotals[keys[2]][mkeys[0]] += d1
                    separatetotals[keys[2]][mkeys[1]] += d1
                    separatetotals[keys[2]][mkeys[2]] += d1

                    separatetotals[mkeys[0]][keys[0]] += d2
                    separatetotals[mkeys[1]][keys[0]] += d2
                    separatetotals[mkeys[2]][keys[0]] += d2
                    separatetotals[mkeys[0]][keys[1]] += d2
                    separatetotals[mkeys[1]][keys[1]] += d2
                    separatetotals[mkeys[2]][keys[1]] += d2
                    separatetotals[mkeys[0]][keys[2]] += d2
                    separatetotals[mkeys[1]][keys[2]] += d2
                    separatetotals[mkeys[2]][keys[2]] += d2

    #  print(all_brackets)
    assert len(all_brackets) == 20
    #  print(all_terms)
    assert len(all_terms) == 10

    for i in sixpoints:
        for j in sixpoints:
            togethertotals[i][j] = round(togethertotals[i][j], 6)
            separatetotals[i][j] = round(separatetotals[i][j], 6)

    for i in sixpoints:
        for j in sixpoints:
            assert togethertotals[i][j] == togethertotals[j][i]
            assert separatetotals[i][j] == separatetotals[j][i]
            assert togethertotals[i][j] == 0.0
            assert separatetotals[i][j] == 0.0
    #  print(togethertotals)
    #  print(separatetotals)

    for i in sixpoints:
        for j in sixpoints:
            if i != j:
                t = septuples[i * 6 + j]
                #  septuples[i * 6 + j] = tuple(round(item, 5) for item in septuples[i * 6 + j])
                #  print(i, j, t, round(sum(t), 5))
                #  assert abs(sum(septuples[i * 6 + j])) < 0.00000001
                t3 = choose(t, 3)
                for s in t3:
                    (a, b, c) = tuple(s)
                    totals = (a + b + c, a + b - c, a - b + c, -a + b + c)
                    totals = tuple([round(item, 8) for item in totals])
                    #  print(totals)
                    for total in totals:
                        if abs(total) < 0.000001:
                            print("!!!!!", s, totals)

    #  30 sets of 4-term GP-relations
    #  30 sets of 6-term GP-relations
    #
#      slabels = []
#      smagnitudes = []
#      p = 2
#      q = 4
#      for i in range(0, 6):
#          for j in range(i + 1, 6):
#              for k in range(j + 1, 6):
#                  #  Uniquely identify the triplet using a key made from prime factors
#                  skeys = (i, j, k)
#                  key = 2 ** (skeys[0] + 1) * 3 ** (skeys[1] + 1) * 5 ** (skeys[2] + 1)
#                  #  Construct the key for the matching face
#                  mkeys = tuple(sorted(list(set(sixpoints).difference([i, j, k]))))
#                  mkey = 2 ** (mkeys[0] + 1) * 3 ** (mkeys[1] + 1) * 5 ** (mkeys[2] + 1)
#                  if not (p in skeys and q in skeys) and not (p in mkeys and q in mkeys):
#                      smagnitudes.append(all_terms[key + mkey] * signpermutation(skeys + mkeys))
#                      slabels.append(label[key + mkey])

#      for i in range(0, 3):
#          for j in range(i + 1, 4):
#              for k in range(j + 1, 5):
#                  for l in range (k + 1, 6):
#                      if smagnitudes[i] + smagnitudes[j] + smagnitudes[k] + smagnitudes[l] == 0:
#                          print(smagnitudes[i], smagnitudes[j], smagnitudes[k], smagnitudes[l])
#                          print(slabels[i], slabels[j], slabels[k], slabels[l])

    print("Calculating number of term-quadruplets that sum to zero.")
    t4 = choose(all_terms.values(), 4)
    print(len(t4))
    quadrupcount = 0
    for q in t4:
        totals = multisign(tuple(q))
        #  print(totals)
        for total in totals:
            if abs(total) < 0.00000001:
                #  print("4-term GP-Relation:", q)
                quadrupcount += 1
    print("4-term GP-Relations=", quadrupcount)

    print("Calculating number of term-sextuplets that sum to zero.")
    t6 = choose(all_terms.values(), 6)
    print(len(t6))
    sextupcount = 0
    for q in t6:
        totals = multisign(tuple(q))
        #  print(totals)
        for total in totals:
            if abs(total) < 0.00000001:
                #  print("6-term GP-Relation:", q)
                sextupcount += 1
    print("6-term GP-Relations=", sextupcount)








    print("Done.")
