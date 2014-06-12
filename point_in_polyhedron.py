"""Code to determine the winding number for a polyhedron around a point, where
that polyhedron is expressed as a (closed, oriented) triangulated surface.

Without loss of generality, we can assume that the point we're interested in is
the origin.

Note: assuming that the given point doesn't lie on any one of the triangles of
the surface, the winding number can be defined as the image of the surface in
the singular homology group H_2 (details needed).

Strategy: we consider the 2-sphere as a CW complex:

 0-cells F and B:  F = (1, 0, 0), B = (-1, 0, 0).
 1-cells L and R:
     L goes from B to F in the x-y plane, through (0, -1, 0)
     R goes from F to B in the x-y plane, through (0, 1, 0)
     Boundary of L is F - B.
     Boundary of R is B - F.
 2-cells U and D:
     U goes through (0, 0, 1).  Boundary is L+R
     D goes through (0, 0, -1). Boundary is -(L+R).
 2-sphere itself is U+D.

Define subsets of R^3 minus the origin:

  P = {(x, y, z) : (x, y, z) > (0, 0, 0) under lexicographic ordering}
  N = {(x, y, z) : (x, y, z) < (0, 0, 0) under lexicographic ordering}

  P intersect (x == 0) is {(y, z) | (y, z) > (0, 0)}, which contains


Map all vertices of the polyhedron to F and B:
  points in P map to F
  points in N map to B

Mapping edges:
  any line segment from a point in N to a point in P must go through
  x == 0 (every point in N has x <= 0, while every point in P has x >= 0).
  There are two cases: (1) both N and P are contained in x == 0, and
  (2) there's a single point of intersection with x == 0.

  In case (1), look at the intersection with y == 0 (the z-axis):
  if it's in P, map to R or -R; if in N, map to L or -L (choose sign
  according to the boundary).

  In case (2), if that single point of intersection lies in P, then
  the edge maps to R or -R (depending on its start and end).  If
  the single point of intersection lies in N, the edge maps to L or -L.

For each triangle:
  map the vertices to F and B.  If all 3 map to the same, nothing to do.
  map the edges to L and R.

"""
import unittest

import numpy


def sign(x):
    """
    Return 1 if x is positive, -1 if it's negative, and 0 if it's zero.

    """
    return (x > 0) - (x < 0)


def ccw2(p1, p2, p3):
    """
    Determine whether the line from p1 to p2 passes to the left or right of p3.

    Return 1 if p1-p2-p3 represents a counterclockwise turn,
    -1 if p1-p2-p3 represents a clockwise turn, and 0 if the three
    points are collinear.

    """
    det = (p2[1] - p3[1]) * (p1[0] - p3[0]) - (p2[0] - p3[0]) * (p1[1] - p3[1])
    return sign(det)


def ccw3(p1, p2, p3, p4):
    """
    Determine whether the counterclockwise triangle p1-p2-p3
    has normal pointing away or towards p4.

    Return 1 if the normal points away from p4, -1 if
    it points towards p4, and 0 if p4 is *contained*
    in the plane of the triangle p1-p2-p3.

    """
    mat = [
        [p1[0] - p4[0], p2[0] - p4[0], p3[0] - p4[0]],
        [p1[1] - p4[1], p2[1] - p4[1], p3[1] - p4[1]],
        [p1[2] - p4[2], p2[2] - p4[2], p3[2] - p4[2]],
    ]
    d = (
        mat[0][0] * mat[1][1] * mat[2][2] +
        mat[0][1] * mat[1][2] * mat[2][0] +
        mat[0][2] * mat[1][0] * mat[2][1] -
        mat[0][0] * mat[1][2] * mat[2][1] -
        mat[0][1] * mat[1][0] * mat[2][2] -
        mat[0][2] * mat[1][1] * mat[2][0])
    return sign(d)


def vertex_chain(v1, origin):
    """
    Map the vertex v1 to the corresponding 0-chain in cellular homology.

    Return a pair (b, f) giving the integer coefficients of B and F
    in C_0 = Z + Z.

    Raise ValueError if v1 == origin.

    """
    d = (sign(v1[0] - origin[0]) or
         sign(v1[1] - origin[1]) or sign(v1[2] - origin[2]))
    if not d:
        raise ValueError("Vertex contains origin")

    if d > 0:
        return numpy.array((0, 1))
    else:
        return numpy.array((1, 0))


def edge_chain(v1, v2, origin):
    """
    Map the edge v1-v2 to the corresponding 1-chain in cellular homology.

    Return a pair (l, r) giving the integer coefficients of L
    and R in C_1 = Z + Z.

    Raise ValueError if the edge goes through the origin.

    """
    edge_boundary = vertex_chain(v2, origin) - vertex_chain(v1, origin)

    assert edge_boundary[0] == -edge_boundary[1]
    if not edge_boundary[0]:
        return numpy.array((0, 0))

    d = ccw2(v1, v2, origin)
    if d == 0:
        if v1[0] != v2[0]:
            d = ccw2((v1[0], v1[2]), (v2[0], v2[2]), (origin[0], origin[2]))
        else:
            d = ccw2((v1[1], v1[2]), (v2[1], v2[2]), (origin[1], origin[2]))
    if not d:
        raise ValueError("Edge contains origin")

    if d * edge_boundary[0] > 0:
        return numpy.array((0, d))
    else:
        return numpy.array((d, 0))


def face_chain(v1, v2, v3, origin):
    """
    Map the triangle v1-v2-v3 to the corresponding 2-chain in cellular
    homology.

    Return result as a pair (d, u) giving coefficients
    of D and U in C_2 = Z + Z.

    """
    face_boundary = sum(
        edge_chain(start, end, origin)
        for start, end in [(v1, v2), (v2, v3), (v3, v1)])

    assert face_boundary[0] == face_boundary[1]
    if not face_boundary[0]:
        return numpy.array((0, 0))

    d = ccw3(v1, v2, v3, origin)
    if not d:
        raise ValueError("Face contains origin")

    if d * face_boundary[0] > 0:
        return numpy.array((0, d))
    else:
        return numpy.array((d, 0))


class Polyhedron(object):
    def __init__(self, triangles, vertex_positions):
        # vertex positions in R^3.
        self.vertex_positions = vertex_positions
        # Indices making up each triangle, counterclockwise
        # around the outside of the face.
        self.triangles = triangles

    def triangle_positions(self):
        """
        Triples of vertex positions.

        """
        for triangle in self.triangles:
            yield tuple(self.vertex_positions[vx] for vx in triangle)

    def winding_number(self, point):
        """Determine whether the given point is inside the polyhedron.

        More precisely, determine the winding number of the polyhedron around
        the point.

        """
        total_face = sum(
            face_chain(v1, v2, v3, point)
            for v1, v2, v3 in self.triangle_positions())
        assert total_face[0] == total_face[1]
        return total_face[0]


# Some sample polyhedra.

tetrahedron = Polyhedron(
    vertex_positions=[
        (-1, -1, -1),
        (-1, 1, 1),
        (1, -1, 1),
        (1, 1, -1),
    ],
    triangles=[
        [0, 1, 3],
        [0, 2, 1],
        [0, 3, 2],
        [1, 2, 3],
    ],
)


octahedron = Polyhedron(
    vertex_positions=[
        (-1, 0, 0),
        (0, -1, 0),
        (0, 0, -1),
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
    ],
    triangles=[
        [0, 2, 1],
        [0, 4, 2],
        [0, 1, 5],
        [0, 5, 4],
        [3, 1, 2],
        [3, 5, 1],
        [3, 2, 4],
        [3, 4, 5],
    ],
)

cube = Polyhedron(
    vertex_positions=[
        (-1, -1, -1),
        (-1, -1, 1),
        (-1, 1, -1),
        (-1, 1, 1),
        (1, -1, -1),
        (1, -1, 1),
        (1, 1, -1),
        (1, 1, 1),
    ],
    triangles=[
        [1, 3, 2],
        [1, 0, 4],
        [1, 5, 7],
        [2, 0, 1],
        [2, 6, 4],
        [2, 3, 7],
        [4, 5, 1],
        [4, 0, 2],
        [4, 6, 7],
        [7, 3, 1],
        [7, 6, 2],
        [7, 5, 4],
    ],
)


class TestPointInPolyhedron(unittest.TestCase):
    def test_tetrahedron(self):
        point = (0, 0, 0)
        poly = tetrahedron
        self.assertEqual(poly.winding_number(point), 1)

        point2 = (1, 1, 1)
        poly = tetrahedron
        self.assertEqual(poly.winding_number(point2), 0)

    def test_cube(self):
        origin = (0, 0, 0)
        self.assertEqual(cube.winding_number(origin), 1)

        # More substantial test, testing a selection of points.
        # Note that we're using floats, but restricting to
        # values that are exactly representable (and for
        # which all the containment calculations are
        # exact).

        def inside(point):
            x, y, z = point
            return -1 < x < 1 and -1 < y < 1 and -1 < z < 1

        def outside(point):
            x, y, z = point
            return x < -1 or x > 1 or y < -1 or y > 1 or z < -1 or z > 1

        # Quarter-integer boundaries from -1.25 to 1.25 inclusive.
        xs = ys = zs = [0.25 * v for v in range(-5, 6)]
        for x in xs:
            for y in ys:
                for z in zs:
                    point = (x, y, z)
                    if inside(point):
                        self.assertEqual(cube.winding_number(point), 1)
                    elif outside(point):
                        self.assertEqual(cube.winding_number(point), 0)
                    else:
                        # Point is on the boundary.
                        with self.assertRaises(ValueError):
                            cube.winding_number(point)

    def test_octahedron(self):
        origin = (0, 0, 0)
        self.assertEqual(octahedron.winding_number(origin), 1)

        # More substantial test, testing a selection of points.
        # Note that we're using floats, but restricting to
        # values that are exactly representable (and for
        # which all the containment calculations are
        # exact).

        def inside(point):
            x, y, z = point
            return abs(x) + abs(y) + abs(z) < 1

        def outside(point):
            x, y, z = point
            return abs(x) + abs(y) + abs(z) > 1

        # Quarter-integer boundaries from -1.25 to 1.25 inclusive.
        xs = ys = zs = [0.25 * v for v in range(-5, 6)]
        for x in xs:
            for y in ys:
                for z in zs:
                    point = (x, y, z)
                    if inside(point):
                        self.assertEqual(octahedron.winding_number(point), 1)
                    elif outside(point):
                        self.assertEqual(octahedron.winding_number(point), 0)
                    else:
                        # Point is on the boundary.
                        with self.assertRaises(ValueError):
                            octahedron.winding_number(point)


if __name__ == '__main__':
    unittest.main()
