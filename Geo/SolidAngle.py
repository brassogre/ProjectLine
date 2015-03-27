'''
Created on Feb 1, 2015

@author: 10034888
'''
import numpy as np

if __name__ == '__main__':

    t = np.random.random((4, 3))
    #  print(t)

#  Examine a very flat tetrahedron   Total Solid angle is 2*pi
#      t[0, 0] = 1.
#      t[0, 1] = 0.0000001
#      t[0, 2] = 0.0000001
#      t[1, 0] = 1.
#      t[1, 1] = 0.01
#      t[1, 2] = 0.01
#      t[2, 0] = 1.
#      t[2, 1] = -0.01
#      t[2, 2] = 0.0000001
#      t[3, 0] = 1.
#      t[3, 1] = 0.0000001
#      t[3, 2] = -0.01

#  Examine a fully symmetric tetrahedron. Total solid angle is 2.30285532191
#      t[0, 0] = -1.0
#      t[0, 1] = -1.0
#      t[0, 2] = +1.0
#      t[1, 0] = -1.0
#      t[1, 1] = +1.0
#      t[1, 2] = -1.0
#      t[2, 0] = +1.0
#      t[2, 1] = -1.0
#      t[2, 2] = -1.0
#      t[3, 0] = +1.0
#      t[3, 1] = +1.0
#      t[3, 2] = +1.0

    indices = (0, 1, 2, 3)
    components = (0, 1, 2)

    #  place the vertices on the unit sphere centered at the origin.
    for i in indices:
        q = 0
        for j in components:
            t[i, j] = t[i, j] * 2.0 - 1.0
            q += t[i, j] ** 2
        distance = np.math.sqrt(q)
        for j in components:
            t[i, j] /= distance

    #  print(t)

    #  Verify that each of the 4 vertices is on the unit sphere.
    for i in indices:
        q = 0
        for j in components:
            q += t[i, j] ** 2
        assert(abs(q - 1.) < 0.00000001)

    #  Compute the three edges (i, k) going from vertex i to vertex k. Each edge has 3 components j
    edges = np.zeros((4, 4, 3))
    dist = np.zeros((4, 4))
    for i in indices:
        for k in set(indices).difference([i]):
            quadrance = 0.0
            for j in components:
                edges[i, k, j] = t[k, j] - t[i, j]
                quadrance += edges[i, k, j] ** 2
            dist[i, k] = np.math.sqrt(quadrance)
            dist[k, i] = dist[i, k]
            assert dist[i, k] <= 2.0
    #  print(edges)

    #  Verify that the edge from i to k is the negative of the edge from k to i
    for i in indices:
        for k in set(indices).difference([i]):
            for j in components:
                assert abs(edges[i, k, j] + edges[k, i, j]) < 0.0000001

    #  For each vertex, calculate the 3 2-d angles on its 3 faces.
    ang = np.zeros((4, 4, 4))
    for i in indices:
        for k1 in set(indices).difference([i]):
            for k2 in set(indices).difference([i, k1]):
                dotproduct = 0.0
                for j in components:
                    dotproduct += edges[i, k1, j] * edges[i, k2, j]
                dotproduct /= (dist[i, k1] * dist[i, k2])
                ang[i, k1, k2] = np.math.acos(dotproduct)



    #  For each vertex, compute the 3 normal vectors associated to its 3 adjoining faces, i.e., 3 pairs of edges.
    normals = np.zeros((4, 4, 3))
    for i in indices:
        for k in set(indices).difference([i]):
            for k1 in set(indices).difference([i, k]):
                for k2 in set(indices).difference([i, k, k1]):
                    #  normal[i, k] is cross product of the edge from i to k1 and from i to k2.
                    #  assert edges[i, k1, 0] and edges[i, k1, 1] and edges[i, k1, 2]
                    #  assert edges[i, k2, 0] and edges[i, k2, 1] and edges[i, k2, 2]
                    normals[i, k, 0] = edges[i, k1, 1] * edges[i, k2, 2] - edges[i, k2, 1] * edges[i, k1, 2]
                    normals[i, k, 1] = -edges[i, k1, 0] * edges[i, k2, 2] + edges[i, k2, 0] * edges[i, k1, 2]
                    normals[i, k, 2] = edges[i, k1, 0] * edges[i, k2, 1] - edges[i, k2, 0] * edges[i, k1, 1]
            #  Adjust the sign of normal so that it points outside the tetrahedron.
            dotproduct = 0.0
            for j in components:
                dotproduct += (t[i, j] - t[k, j]) * normals[i, k, j]
            if dotproduct > 0.0:
                #  Flip the sign of the normal.
                for j in components:
                    normals[i, k, j] *= -1.0
    #  print(normals)

    #  Make the normal vectors unit length
    for i in indices:
        for k in set(indices).difference([i]):
            q = 0 #  set the initial quadrance at 0
            for j in components:
                assert normals[i, k, j]
                q += normals[i, k, j] ** 2
            distance = np.math.sqrt(q)
            for j in components:
                normals[i, k, j] /= distance
    #  print(normals)

    #  Verify that each normal vector has unit length
    for i in indices:
        for k in indices:
            if i == k:
                assert normals[i, k, 0] == 0
                assert normals[i, k, 1] == 0
                assert normals[i, k, 2] == 0
            else:
                assert abs(normals[i, k, 0] ** 2 + normals[i, k, 1] ** 2 + normals[i, k, 2] ** 2 - 1.0) < 0.0000001

    #  Verify that each normal points out of the tetrahedron
    for i in indices:
        for k in set(indices).difference([i]):
            dotproduct = 0.0
            for j in components:
                dotproduct += (t[i, j] - t[k, j]) * normals[i, k, j]
            if dotproduct > 0.0:
                print("Need to flip sign for ", end='')
                print(i, k)

    #  Calculate the angles between the normals
    angles = np.zeros((4, 4))
    for i in indices:
        for k in set(indices).difference([i]):
            innerindices = set(indices).difference([i, k])
            k2 = innerindices.pop()
            k1 = innerindices.pop()
            assert len(innerindices) == 0
            assert k1 != k2
            assert k1 != k and k2 != k
            assert k1 != i and k2 != i
            dotproduct = 0
            for j in components:
                assert normals[i, k1, j]
                assert normals[i, k2, j]
                dotproduct += normals[i, k1, j] * normals[i, k2, j]
            angles[i, k] = np.math.acos(-dotproduct)
            assert angles[i, k] > 0

    #  Verify that the matrix of angles is symmetric
    for i in indices:
        for k in indices:
            assert abs(angles[i, k] - angles[k, i]) < 0.00000001

    #  Calculate the standard solid angles.
    solid = np.full((4), -np.math.pi)
    sol = [0, 0, 0, 0]
    for i in indices:
        assert solid[i] == -np.math.pi
        for k in set(indices).difference([i]):
            solid[i] += angles[i, k]
#          for k1 in set(indices).difference([i]):
#              for k2 in set(indices).difference([i, k1]):
#                  key = set(indices).difference([i, k1, k2]).pop() + k1 + k2 - 3
#                  sol[key] = solid[i]
    print(solid)
    #  print(sol)

    #  Sum up all 4 solid angles
    total = 0
    for i in indices:
        total += solid[i]
    #  print(total)

    #  Sum twice the 6 dihedral angles minus 4pi.
    othertotal = -4 * np.math.pi
    for i in indices:
        for k in range(i + 1, 4):
            othertotal += 2 * angles[i, k]
    #  print(othertotal)
    assert abs(total - othertotal) < 0.00000001

    #  Investigate possible generalizations of curvature (like Spread/Quadrance)
#      curvature2d = [(dist[2, 3], angles[0, 1]),
#        (dist[1, 3], angles[0, 2]),
#        (dist[1, 2], angles[0, 3]),
#        (dist[0, 3], angles[1, 2]),
#        (dist[0, 2], angles[1, 3]),
#        (dist[0, 1], angles[2, 3])]
#
#      curvature3 = [(dist[2, 3], angles[0, 1], angles[0, 1]),
#        (dist[1, 3], dist[0, 2], angles[0, 2]),
#        (dist[1, 2], dist[0, 3], angles[0, 3]),
#        (dist[0, 3], dist[1, 2], angles[1, 2]),
#        (dist[0, 2], dist[1, 3], angles[1, 3]),
#        (dist[0, 1], dist[2, 3], angles[2, 3])]

#      print([ np.math.sin(a) / x for (x, a) in curvature2d ])
#      print([ np.math.sin(a) * y / x for (x, y, a) in curvature3 ])
#      print([ np.math.sin(a) * y / x ** 2 for (x, y, a) in curvature3 ])
#      print([ np.math.sin(a) ** 2 * y / x for (x, y, a) in curvature3 ])
#      print([ np.math.sin(a) ** 2 * y / x ** 2 for (x, y, a) in curvature3 ])

    #  Create a list of the curvatures (sin(x)/x values) for each face
    dictmap = { }
    for i in indices:
        for k1 in set(indices).difference([i]):
            for k2 in set(indices).difference([i, k1]):
                key = set(indices).difference([i, k1, k2]).pop()
                val = np.math.sin(ang[i, k1, k2]) / dist[k1, k2]
                lll = dictmap.get(key, [])
                lll.append(val)
                dictmap[key] = lll
    curvs = []
    for item in dictmap:
        lll = dictmap[item]
        firstelem = lll[0]
        curvs.append(lll[0])
        for elem in lll:
            assert abs(firstelem - elem) < 0.0000001
    print(curvs)

    #  Create a list of the circumcircle areas for each face
    circumareas = list(map(lambda x : np.math.pi / (4 * x), curvs))
    print(circumareas)

    #  Create a list of the triangle areas for each face
    dictmap = { }
    triareas = [0, 0, 0, 0]
    for i in indices:
        for k1 in set(indices).difference([i]):
            for k2 in set(indices).difference([i, k1]):
                key = set(indices).difference([i, k1, k2]).pop()
                val = np.math.sin(ang[i, k1, k2]) * dist[i, k1] * dist[i, k2] / 2.0
                #  print(key, val)
                tr = dictmap.get(key, None)
                if tr:
                    assert abs(tr - val) < 0.0000001
                else:
                    dictmap[key] = val
    for item in dictmap:
        triareas[item] = dictmap[item]
    print(triareas)
    print()

    prodcurv = [0, 0, 0, 0]
    for i in indices:
        for k1 in set(indices).difference([i]):
            for k2 in set(indices).difference([i, k1]):
                key = set(indices).difference([i, k1, k2]).pop()
                prodcurv[key] = curvs[i] * curvs[k1] * curvs[k2]

#      funcs = { np.math.sin, np.math.log }
#      for f in funcs:
#          print(list(map(f, triareas)))

    #  Calculate the circumquadrances

    #  Calculate some cross products
