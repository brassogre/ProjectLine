'''
Created on Mar 22, 2015

@author: Kutach
'''
import fractions

def reduce(seqList):
    """
    Eliminates any factors that are common to all the numbers making up this proportion.
    This is used primarily at the last step of initializing Ratios.
    """
    n = len(seqList)
    if n > 0:
        gcd = seqList[-1]
        #  Search through all the terms except the last (the denominator) for the
        #  greatest common divisor
        for v in seqList:
            gcd = fractions.gcd(gcd, v)

        #  Now divide out the greatest common divisor to keep the numbers small.
        if gcd > 1:
            for i in range(0, n):
                assert seqList[i] % gcd == 0
                seqList[i] //= gcd
    return seqList

class Ratio:
    def __init__(self, *args):
        if len(args) == 0:
            raise ValueError() #  Not allowed to create a Ratio of undefined value.
        if len(args) == 1 and isinstance(args[0], Ratio):
            values = args[0].values
        else:
            values = args
        newvalues = []
        for i in values:
            newvalues.append(int(i))
        if len(values) == 1:
            newvalues.append(1)
        self.values = tuple(reduce(newvalues)) #  Make this ratio immutable

#      def __int__(self):
#          if self.values[-1] != 1 or len(self.values) > 2:
#              raise Exception()
#          else:
#              return self.values[0]

    def __repr__(self):
        s = "<"
        for i in self.values:
            s += str(i) + ":"
        s = s[:-1] + ">" #  Eliminate final colon.
        return s

    def multiply(self, n):
        denominator = self.values[-1]
        newvalue = []
        gcd = fractions.gcd(n, denominator)
        if gcd > 1:
            denominator //= gcd
            n //= gcd

        #  Multiply all the terms except the last (the denominator) by n
        for v in self.values[:-1]:
            newvalue.append(v * n)
        newvalue.append(denominator)

        #  Return the new value
        return Ratio(*newvalue)

    def divide(self, n):
        newvalue = list(self.values)
        if isinstance(n, Ratio): #  TODO: Figure out how to divide by multi-valued Ratios...
            newvalue[-1] = n.multiply(newvalue[-1])
        else:
            newvalue[-1] = newvalue[-1] * n
        #  Return the new value
        return Ratio(*newvalue)

    def samesize(self, r):
        return isinstance(r, Ratio) and len(r.values) == len(self.values)

    def add(self, r):
        if not self.samesize(r):
            raise ValueError()
        newvalues = []
        myden = self.values[-1]
        otherden = r.values[-1]
        for i in range(0, len(self.values) - 1):
            newvalues.append(self.values[i] * otherden + r.values[i] * myden)

        newvalues.append(myden * otherden)
        #  Return the new value
        return Ratio(*newvalues)

    def subtract(self, r):
        if not self.samesize(r):
            raise ValueError()
        newvalues = []
        myden = self.values[-1]
        otherden = r.values[-1]
        for i in range(0, len(self.values) - 1):
            newvalues.append(self.values[i] * otherden - r.values[i] * myden)

        newvalues.append(myden * otherden)
        #  Return the new value
        return Ratio(*newvalues)

    def len(self):
        return len(self.values)

def zeros(n):
    matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(0)
        matrix.append(row)
    return matrix

def ones(n):
    matrix = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(1)
        matrix.append(row)
    return matrix

def simplistic_determinant(A):
    """
    This implements a naive slow algorithm just as a double check that the Crout-based algorithm works.
    """
    n = len(A)
    if n == 0:
        return None
    #  Base case
    if n == 1:
        return A[0][0]

    #  Inductive case
    det = 0
    for row in range(0, n):
        coefficient = A[row][0]
        sign = (-1) ** row
        minor = zeros(n - 1)
        rshift = 0
        for r in range(0, n - 1):
            if r == row:
                rshift = 1
            for c in range(0, n - 1):
                minor[r][c] = A[r + rshift][c + 1]
        det += sign * coefficient * simplistic_determinant(minor)
#      if n == 2:
#          print(A, det)
    return det

def determinant(A):
    """
    Returns the detrminant of the matrix A, whose elements are assumed to be integers.
    The computation uses Crout's Algorithm to perform LU decomposition on A.

    INPUT:
        A: a sequence of sequences representing a square matrix

    OUTPUT:
        an int representing the determinant of A
    """
    #  This is Crout's Algorithm.
    #  U will remain zero in the lower left entries, and will be 1's along the diagonal.
    #  L will remain zero in the upper right entries.
    #  A = U L
    n = len(A)
    L = zeros(n) #  Initialize with zeros  Numerators for the lower triangular matrix
    U = zeros(n) #  Initialize with zeros  Numerators for the upper triangular matrix
    DL = ones(n) #  Initialize with zeros  Denominators for the lower triangular matrix
    DU = ones(n) #  Initialize with zeros  Denominators for the upper triangular matrix
    #  L = [[0] * n] * n #  Does not work because it initializes the matrix with references to the same lists
    #  U = [[0] * n] * n #  Does not work because it initializes the matrix with references to the same lists
    for j in range(0, n):
        assert len(A[j]) == n
        U[j][j] = 1             #  set the diagonal entries of U to 1
        for i in range(j, n):  #  starting at L[j][j], solve j-th column of L
            tempL = A[i][j]
            tempDL = 1 #  Temporary denominator for the lower triangular matrix
            for k in range(0, j):
                assert DL[i][k] != 0
                assert DU[k][j] != 0
                tempL = tempL * DL[i][k] * DU[k][j] - tempDL * L[i][k] * U[k][j]
                tempDL = tempDL * DL[i][k] * DU[k][j]
            L[i][j] = tempL
            DL[i][j] = tempDL
        for i in range(j + 1, n):#  starting at U[j][j+1], solve j-th row of U
            tempU = A[j][i]
            tempDU = 1 #  Temporary denominator for the upper triangular matrix
            for k in range(0, j):
                assert DU[k][i] != 0
                assert DL[j][k] != 0
                tempU = tempU * DU[k][i] * DL[j][k] - tempDU * L[j][k] * U[k][i]
                tempDU = tempDU * DU[k][i] * DL[j][k]
            U[j][i] = tempU * DL[j][j]
            if L[j][j] == 0:
                assert simplistic_determinant(A) == 0
                return 0 #  The determinant is zero, so avoid dividing by zero by short circuiting the computation.
            DU[j][i] = tempDU * L[j][j]

    #  Now calculate the determinant by multiplying the diagonal entries of the lower-left triangular matrix
    num = 1
    den = 1
    for i in range(0, n):
        assert U[i][i] == 1
        for j in range(0, i):
            assert U[i][j] == 0
        for j in range(i + 1, 3):
            assert L[i][j] == 0
        num *= L[i][i]
        den *= DL[i][i]
    #  Now divide the denominator den from the numerator.
    #  The numerator should evenly divide (assuming the input matrix A only had ints),
    #  Was having trouble with den equaling 0. Fixed it by returning zero if any L[j][j] was 0. See above.
    assert den != 0
    det1 = num // den
    det2 = simplistic_determinant(A)
    if det1 != det2:
        print("Mismatch! ", det1, det2)
    return det2

class ScalarVariable:
    """
    Represents scalar variables using a label and subscript
    """
    def __init__(self, label, subscript=None, *args, **keyargs):
        if not isinstance(label, str):
            raise TypeError()
        self.label = label[0]
        self.subscript = subscript

    def __repr__(self):
        value = "$" + self.label
        if self.subscript != None:
            value += "_{" + str(self.subscript) + "}"
        value += "$"
        return value

class PointVariable(ScalarVariable):
    """
    Represents projective points using coordinates
    """
    def __init__(self, label, subscript=None, dimension=None, coordinates=None, *args, **keyargs):
        super().__init__(label, subscript)
        self.coordinates = None
        if coordinates:
            if isinstance(coordinates, Ratio):
                copy_of_coordinates = list(coordinates.values)
            else:
                copy_of_coordinates = list(coordinates)
#              if not isinstance(coordinates, collections.abc.Sequence):
#                  raise TypeError() #  coordinates need to be a Ratio or a Sequence
            d = len(copy_of_coordinates)
            if dimension:
                if d != dimension:
                    raise TypeError() #  Explicit dimension needs to match number of coordinates passed in.
            else:
                dimension = d
            self.coordinates = tuple(copy_of_coordinates)
        if not dimension or not isinstance(dimension, int):
            raise TypeError() #  Construction requires either explicit Integral dimension or coordinates

        self.dimension = dimension

    def __repr__(self):
        value = super().__repr__()
        if self.coordinates:
            value = value[:-1] + "=" + str(self.coordinates) + "$"
        return value

    def has_coordinates(self):
        return self.coordinates != None

class BracketVariable:
    """
    Represents projective bracket monomials using PointVariables as entries.
    The BracketVariable keeps track of the anti-commutativity of the PointVariables in the BracketVariable.
    """
    def __init__(self, *args, **keyargs):
        dim = None #  Let dimension be initially unknown so that it can be set by the first PointVariable passed.
        entries = [] #  PointVariables will be collected in this list, then later saved as a tuple.
        b_can_calculate_value = True #  boolean flag to set equal to false if any of the PointVariables does not have coordinates.
        for item in args:
            if not isinstance(item, PointVariable):
                raise TypeError() #  Brackets need to be composed only of PointVariables.
            if dim:
                if item.dimension != dim:
                    raise TypeError() #  All passed PointVariables need to be of the same dimension.
            else:
                dim = item.dimension #  Set dim to the dimension of the first PointVariable seen.
            entries.append(item)
            b_can_calculate_value &= item.has_coordinates()
        self.entries = tuple(entries)

        self.value = None
        if b_can_calculate_value:
#              print(self.entries)
#              for pv in self.entries:
#                  assert isinstance(pv, PointVariable)
            M = tuple([pv.coordinates for pv in self.entries])
            #  print(M)
            self.value = determinant(M)

    def __repr__(self):
        value = "["
        for pv in self.entries:
            value += str(pv)
            if pv != self.entries[-1]:
                value += ", "
        value += "]"
        if self.value != None:
            value += "=" + str(self.value) + " "
        return value

#  TODO: Create BracketMonomial

#  TODO: Create BracketPolynomial

#  TODO: Check whether a rational parameterization of the sphere will serve as better than the standard embedding of the projective plane.
#  Rational parameterization of a circle in homogeneous coordinates is (2t, 1-t^2, 1+t^2)
#  Rational parameterization of a sphere in homogeneous coordinates is (2s, 2t, 1-s^2-t^2, 1+s^2+t^2)


if __name__ == '__main__':
#      a = ScalarVariable("ab", 2)
#      b = ScalarVariable("b", 0)
#      C = ScalarVariable("C")
#      print(a, b, C)

#      s = PointVariable("s", coordinates=(1, 2, 3))
#      t = PointVariable("t", coordinates=(1, 2))
#      r = PointVariable("r", dimension=3)
#      p = PointVariable("p", 6, coordinates=(1, 2, 3))
#      q = PointVariable("q", dimension=2, coordinates=(6, 7))
#    print(s, t, r, p, q)

    a = PointVariable("a", coordinates=(-1, 3, 1))
    b = PointVariable("b", coordinates=(1, 1, 1))
    c = PointVariable("c", coordinates=(2, 2, 1))
    d = PointVariable("d", coordinates=(101, 101, 1))
    br1 = BracketVariable(a, b, c)
    br4 = BracketVariable(c, b, d)
    br2 = BracketVariable(a, b, d)
    br3 = BracketVariable(a, d, c)
    print(br1.value + br4.value, br2.value + br3.value)
    print(br1.value, br4.value * br2.value * br3.value)

#      test = {0:[[2, 0, 0], [0, 3, 0], [0, 0, 5]], \
#              1:[[1, 2, 3], [2, 3, 4], [5, 6, 7]], \
#              2:[[-3, 0, 0], [0, -3, 0], [0, 0, -3]], \
#              3:[[2, 3, 5], [7, 11, 13], [17, 19, 23]], \
#              4:[[2, 3, 5], [0, 0, 0], [17, 19, 23]], \
#              5:[[2, 3, 5], [7, -19, 0], [-7, 19, 0]] }
#
#      for key in test:
#          det = determinant(test[key])
#          print("Determinant = ", det)

