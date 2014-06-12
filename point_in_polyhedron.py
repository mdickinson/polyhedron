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


def sign(x):
    """
    Return 1 if x is positive, -1 if it's negative, and 0 if it's zero.

    """
    return (x > 0) - (x < 0)


def ccw3(p1, p2, p3, p4):
    """
    Determine whether the counterclockwise triangle p1-p2-p3
    has normal pointing away or towards p4.

    Return 1 if the normal points away from p4, -1 if
    it points towards p4, and 0 if p4 is *contained*
    in the plane of the triangle p1-p2-p3.

    """
    mat = [
        [p1[0] - p4[0], p2[0] - p4[0], p3[0]- p4[0]],
        [p1[1] - p4[1], p2[1] - p4[1], p3[1]- p4[1]],
        [p1[2] - p4[2], p2[2] - p4[2], p3[2]- p4[2]],
    ]
    d = (
        mat[0][0] * mat[1][1] * mat[2][2] +
        mat[0][1] * mat[1][2] * mat[2][0] +
        mat[0][2] * mat[1][0] * mat[2][1] -
        mat[0][0] * mat[1][2] * mat[2][1] -
        mat[0][1] * mat[1][0] * mat[2][2] -
        mat[0][2] * mat[1][1] * mat[2][0])
    return sign(d)


def ccw(p1, p2, p3):
    """
    Return 1 if p1-p2-p3 represents a counterclockwise turn,
    -1 if p1-p2-p3 represents a clockwise turn, and 0 if the three
    points are collinear.

    """
    det = (p2[1] - p3[1]) * (p1[0] - p3[0]) - (p2[0] - p3[0]) * (p1[1] - p3[1])
    if det == 0:
        return 0
    elif det > 0:
        return 1
    else:
        return -1


def classify_vertex(v1, origin):
    """
    Classify vertex with respect to the origin.

    Returns (1, 0) for points in N, (0, 1) for points in P.

    Raise ValueError if v1 == origin.

    """
    if v1 == origin:
        raise ValueError("cannot classify origin")
    return (1, 0) if v1 < origin else (0, 1)


def classify_edge(v1, v2, origin):
    """
    Classify edge from point v1 in N to point v2 in P,
    relative to the origin.

    Return a pair (l, r) giving the integer coefficients of L
    and R in the result.

    """
    if v2 < origin < v1:
        x, y = classify_edge(v2, v1, origin)
        return -x, -y

    if not v1 < origin < v2:
        raise ValueError("edge not between N and P")

    # CCW rule ignoring 3rd coordinate.
    y = (v2[0] - origin[0]) * (v1[1] - origin[1]) - (v2[1] - origin[1]) * (v1[0] - origin[0]) 
    if y > 0:
        return (0, -1)
    elif y < 0:
        return (1, 0)

    # Case where v1, v2 and origin are collinear in x-y plane, so
    # Find z-coordinate of point of intersection.

    # Note that it's not possible for v1 and v2 to both lie on the z-axis,
    # since then the edge would pass through the origin.

    # Let's assume x-coordinates are distinct.
    if v1[0] < v2[0]:
        z = (v2[0] - origin[0]) * (v1[2] - origin[2]) - (v2[2] - origin[2]) * (v1[0] - origin[0])
    else:
        z = (v2[1] - origin[1]) * (v1[2] - origin[2]) - (v2[2] - origin[2]) * (v1[1] - origin[1])

    if z > 0:
        return (0, -1)
    elif z < 0:
        return (1, 0)
    else:
        raise ValueError("edge through origin")


def classify_triangle(v1, v2, v3, origin):
    """Classify triangle, assuming that v1-v2-v3-v1 gives a cycle in the 1-complex.

    Assume that the cycle v1-v2-v3-v1 maps to L+R, which implies in turn that
    v1, v2 and v3 are non-collinear, and that the point of intersection exists
    and is unique.

    """
    # Figure out whether the plane containing v1, v2 and v3
    # lies above the origin (in the sense of the z-axis) or not.

    # Equation of plane: (r - origin) . n = d, where:
    #
    #   normal n = (v2 - v1) x (v3 - v1)
    #          d = (v2 - v1) x (v3 - v1) . (v1 - origin)
    #
    # We only need the z-component of the normal.
    #
    # Now find t such that r = origin + t * k satisfies the equation:
    #
    #    (t * k) . n = d
    #
    # that is, 
    #
    #    t * nz = d
    #
    # Thus the sign of t matches that of nz * d, provided that nz * d is nonzero.
    #
    # d is zero iff the face contains the origin; in this case, raise
    # an exception.
    #
    # If nz == 0 then the triangle is vertical.  This means that the entire
    # edge cycle projects down onto a line-segment.  The only way that this line can
    # be a nontrivial cycle is if the face contains the origin, and in that
    # case d will also be zero.

    # XXX Computing nz is redundant: we already have the classification
    # of the edge cycle.

    nz = (v2[0] - v1[0]) * (v3[1] - v1[1]) - (v2[1] - v1[1]) * (v3[0] - v1[0])
    nz = sign(nz)
    # d is the determinant of the 3x3 matrix |v1-origin v2-origin v3-origin|,
    # or equivalently, of the 4x4 matrix
    #
    #   |v1 v2 v3 origin|
    #   | 1  1  1   1   |

    d = ccw3(v1, v2, v3, origin)
    if not d:
        raise ValueError("Face contains origin")

    assert nz

    t_sign = d * nz
    if t_sign == 0:
        raise ValueError("Face contains origin")
    elif t_sign > 0:
        return (0, 1)
    else:
        return (-1, 0)


class Polyhedron(object):
    def __init__(self, triangles, vertex_positions):
        # vertex positions in R^3.
        self.vertex_positions = vertex_positions
        # Indices making up each triangle, counterclockwise
        # around the outside of the face.
        self.triangles = triangles

    def winding_number(self, point):
        """Determine whether the given point is inside the polyhedron.

        More precisely, determine the winding number of the polyhedron around
        the point.

        """
        face_contributions = []
        for triangle in self.triangles:
            # Classify positions.  If all 3 lie in the same half-space,
            # there's nothing of interest.
            v1, v2, v3 = [self.vertex_positions[v] for v in triangle]
            p1, p2, p3 = [classify_vertex(self.vertex_positions[v], point) for v in triangle]

            # Now classify edges from N to P (or P to N).
            edge_classifications = []

            pairs = [(v1, v2), (v2, v3), (v3, v1)]
            cpairs = [(p1, p2), (p2, p3), (p3, p1)]
            for pair, cpair in zip(pairs, cpairs):
                c1, c2 = cpair
                if c1 == c2:
                    continue
                vx1, vx2 = pair
                ec = classify_edge(vx1, vx2, point)
                edge_classifications.append(ec)

            total_edge = sum(e[0] for e in edge_classifications), sum(e[1] for e in edge_classifications)
            # Those with total_edge != (0, 0) contribute to the winding number;
            # those with total_edge == (0, 0) do not.
            if total_edge[0]:
                # Now classify face: we want to know whether the plane
                # through v1, v2 and v3 intersects the z-axis above or below the origin.
                d, u = classify_triangle(v1, v2, v3, point)
                d, u = total_edge[0] * d, total_edge[0] * u
                face_contributions.append((d, u))

        total_face = sum(f[0] for f in face_contributions), sum(f[1] for f in face_contributions)
        assert total_face[0] == total_face[1]
        return total_face[0]


# Some sample polyhedra.

tetrahedron = Polyhedron(
    vertex_positions = [
        (-1, -1, -1),
        (-1, 1, 1),
        (1, -1, 1),
        (1, 1, -1),
    ],
    triangles = [
        [0, 1, 3],
        [0, 2, 1],
        [0, 3, 2],
        [1, 2, 3],
    ],
)


octahedron = Polyhedron(
    vertex_positions = [
        (-1, 0, 0),
        (0, -1, 0),
        (0, 0, -1),
        (1, 0, 0),
        (0, 1, 0),
        (0, 0, 1),
    ],
    triangles = [
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
    vertex_positions = [
        (-1, -1, -1),
        (-1, -1, 1),
        (-1, 1, -1),
        (-1, 1, 1),
        (1, -1, -1),
        (1, -1, 1),
        (1, 1, -1),
        (1, 1, 1),
    ],
    triangles = [
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

    def test_classify_edge(self):
        # classify edges of a cube.
        origin = (0, 0, 0)
        # vertices: v1, v3, v5, v7 lie in N, the others in P.
        v1 = (-1, -1, -1)
        v2 = (1, -1, -1)
        v3 = (-1, 1, -1)
        v4 = (1, 1, -1)
        v5 = (-1, -1, 1)
        v6 = (1, -1, 1)
        v7 = (-1, 1, 1)
        v8 = (1, 1, 1)

        # Some midpoints.
        m12 = (0, -1, -1)  # in N
        m34 = (0, 1, -1)  # in P
        m56 = (0, -1, 1)  # in N
        m78 = (0, 1, 1)  # in P

        # Edges from N to P.
        self.assertEqual(classify_edge(v1, v2, origin), (1, 0))
        self.assertEqual(classify_edge(v3, v4, origin), (0, -1))
        self.assertEqual(classify_edge(v5, v6, origin), (1, 0))
        self.assertEqual(classify_edge(v7, v8, origin), (0, -1))
        # Edges from P to N.
        self.assertEqual(classify_edge(v2, v1, origin), (-1, 0))
        self.assertEqual(classify_edge(v4, v3, origin), (0, 1))
        self.assertEqual(classify_edge(v6, v5, origin), (-1, 0))
        self.assertEqual(classify_edge(v8, v7, origin), (0, 1))

        # Edges from N to N or P to P should raise.
        vn = [v1, v3, v5, v7]
        vp = [v2, v4, v6, v8]
        for start in vn:
            for end in vn:
                with self.assertRaises(ValueError):
                    classify_edge(start, end, origin)
        for start in vp:
            for end in vp:
                with self.assertRaises(ValueError):
                    classify_edge(start, end, origin)

        # Some face diagonals.
        # Line segment from v1 (in N) to v4 (in P) should map
        # to either (1, 0) or (0, -1).  The projection to the x-y plane
        # goes through the origin, so the first CCW test is inconclusive.
        self.assertEqual(classify_edge(v1, v4, origin), (1, 0))
        self.assertEqual(classify_edge(v4, v1, origin), (-1, 0))

        # Line segment from v3 (in N) to v2 (in P) should do the same.
        self.assertEqual(classify_edge(v3, v2, origin), (1, 0))
        self.assertEqual(classify_edge(v2, v3, origin), (-1, 0))

        # and same through the positive part of z-axis.
        self.assertEqual(classify_edge(v5, v8, origin), (0, -1))
        self.assertEqual(classify_edge(v8, v5, origin), (0, 1))
        self.assertEqual(classify_edge(v7, v6, origin), (0, -1))
        self.assertEqual(classify_edge(v6, v7, origin), (0, 1))

        # Troublesome edges, going through the origin after projection
        # to the x-y plane.  Classification depends on the point
        # at which the edge intersects the z-axis.
        self.assertEqual(classify_edge(m12, m34, origin), (1, 0))
        self.assertEqual(classify_edge(m34, m12, origin), (-1, 0))
        self.assertEqual(classify_edge(m56, m78, origin), (0, -1))
        self.assertEqual(classify_edge(m78, m56, origin), (0, 1))

        # Body diagonals should raise
        with self.assertRaises(ValueError):
            classify_edge(v1, v8, origin)
        with self.assertRaises(ValueError):
            classify_edge(v2, v7, origin)
        with self.assertRaises(ValueError):
            classify_edge(v3, v6, origin)
        with self.assertRaises(ValueError):
            classify_edge(v4, v5, origin)


if __name__ == '__main__':
    unittest.main()
