"""Simple point-in-polygon algorithm based on winding number, with robustness
depending only on the underlying arithmetic.

We've got a closed, possibly non-simple polygon described as a list of vertices
in R^2, and we're given a point that doesn't lie directly on the path of the
polygon.  We'd like to compute the winding number of the polygon around the
point.

There are two sources of difficulty: (1) dealing with numerical errors that
might result in an incorrect answer, and (2) dealing with degenerate cases.
We'll ignore the numerical issues for the moment.

Strategy: without loss of generality, let's place the point at the origin.
Divide the remainder of the plane (i.e., R^2 minus the origin) into two
halves, L and R, defined as follows:

    L = {(x, y) | x < 0 or x == 0 and y < 0}

    R = {(x, y) | x > 0 or x == 0 and y > 0}

That is, R contains all points with argument in the half-closed interval
(-pi/2, pi/2], and L contains all others.  Note that with these definitions, L
and R are both convex: a line segment between two points in R lies entirely in
R, and similarly for L.  In particular, a line segment between two points can
only pass through the origin if one of those points is in L and the other in R.

Also, with these conventions it's simple to test whether a point is in L: as a
Python tuple, (x, y) is in L iff (x, y) < (0, 0).

Now the idea is that we follow the edges of the polygon, keeping track of how
many times we move between L and R.  For each move from L to R (or vice versa),
we also need to compute whether the edge passes *above* or *below* the origin,
to compute its contribution to the total winding number.  From the comment
above, we can safely ignore all edges that lie entirely within either L or R.

"""

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


def half_turn(point1, point2, origin):
    """
    Return the contribution to the total winding number about 'origin' from a
    single edge, from point1 to point2.

    """
    if point1 == origin or point2 == origin:
        raise ValueError("Vertex at origin")

    # If both points lie in the same half-plane relative to the
    # origin, return 0.
    if (point1 < origin) == (point2 < origin):
        return 0

    turn = ccw(point1, point2, origin)
    if not turn:
        raise ValueError("Segment goes through origin")
    return turn


def winding_number(polygon, origin):
    """
    Compute the (counterclockwise) winding number of a polygon around a point.

    Raise ValueError if the point lies directly on the path
    of the polygon.

    """
    edges = zip(polygon, polygon[1:] + polygon[:1])
    return sum(
        half_turn(point1, point2, origin) for point1, point2 in edges) // 2


import unittest

class TestWindingNumber(unittest.TestCase):
    def test_simple_square(self):
        square = [
            (1.0, -1.0),
            (1.0, 1.0),
            (-1.0, 1.0),
            (-1.0, -1.0),
        ]
        origin = (0.0, 0.0)
        self.assertEqual(winding_number(square, origin), 1)

    def test_double_square(self):
        square = [
            (1.0, -1.0),
            (1.0, 1.0),
            (-1.0, 1.0),
            (-1.0, -1.0),
        ] * 2
        origin = (0.0, 0.0)
        self.assertEqual(winding_number(square, origin), 2)

    def test_clockwise_square(self):
        square = [
            (1.0, -1.0),
            (1.0, 1.0),
            (-1.0, 1.0),
            (-1.0, -1.0),
        ][::-1]
        origin = (0.0, 0.0)
        self.assertEqual(winding_number(square, origin), -1)

    def test_various_points_in_square(self):
        square = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        test_points = []
        for x in range(-3, 4):
            for y in range(-3, 4):
                test_points.append((0.5*x, 0.5*y))

        for point in test_points:
            x, y = point
            if -1 < x < 1 and -1 < y < 1:
                # Point is inside.
                self.assertEqual(winding_number(square, point), 1)
            elif x < -1 or x > 1 or y < -1 or y > 1:
                # Point outside.
                self.assertEqual(winding_number(square, point), 0)
            else:
                with self.assertRaises(ValueError):
                    winding_number(square, point)

    def test_aitch(self):
        aitch = [
            (0, 0),
            (1, 0),
            (1, 1),
            (2, 1),
            (2, 0),
            (3, 0),
            (3, 3),
            (2, 3),
            (2, 2),
            (1, 2),
            (1, 3),
            (0, 3),
        ]

        test_points = [
            (0.5*x, 0.5*y)
            for y in range(-1, 8)
            for x in range(-1, 8)
        ]

        # * for boundary, '.' for outside, 'o' for inside.
        template = """\
.........
.***.***.
.*o*.*o*.
.*o***o*.
.*ooooo*.
.*o***o*.
.*o*.*o*.
.***.***.
.........
"""
        expected = ''.join(template.strip().split())

        assert len(expected) == len(test_points)
        for point, point_type in zip(test_points, expected):
            if point_type == '.':
                self.assertEqual(winding_number(aitch, point), 0)
            elif point_type == 'o':
                self.assertEqual(winding_number(aitch, point), 1)
            else:
                with self.assertRaises(ValueError):
                    winding_number(aitch, point)


def test_main(exit=True):
    unittest.main(exit=exit)


def show_winding_number():
    import random

    import numpy
    from matplotlib import pyplot

    vertex_count = 7
    polygon = [
        (random.random(), random.random())
        for _ in range(vertex_count)
    ]

    # Select a collection of points, and compute their winding numbers.
    xs = numpy.random.rand(10000)
    ys = numpy.random.rand(10000)
    ws = []
    for x, y in zip(xs, ys):
        point = (x, y)
        ws.append(winding_number(polygon, point))
    ws = numpy.array(ws, dtype=float)

    xxs, yys = zip(*polygon)
    pyplot.plot(xxs + (xxs[0],), yys + (yys[0],))
    sc = pyplot.scatter(xs, ys, c=ws)
    pyplot.colorbar(sc)
    pyplot.show()




if __name__ == "__main__":
    test_main(exit=False)
    show_winding_number()
