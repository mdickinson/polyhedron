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

from polyhedron import Polyhedron


# Some sample polyhedra.

# Regular tetrahedron.
tetrahedron = Polyhedron(
    vertex_positions=[
        (0, 0, 0),
        (0, 1, 1),
        (1, 0, 1),
        (1, 1, 0),
    ],
    triangles=[
        [0, 1, 3],
        [0, 2, 1],
        [0, 3, 2],
        [1, 2, 3],
    ],
)


def tetrahedron_classify(point):
    x, y, z = point
    if x + y + z < 2 and x + y > z and x + z > y and y + z > x:
        return "inside"
    if x + y + z <= 2 and x + y >= z and x + z >= y and y + z >= x:
        return "boundary"
    return "outside"


# Regular octahedron, with vertices on the axes.
octahedron = Polyhedron(
    vertex_positions=[
        (-1, 0, 0), (0, -1, 0), (0, 0, -1),
        (1, 0, 0), (0, 1, 0), (0, 0, 1),
    ],
    triangles=[
        [0, 2, 1], [0, 4, 2], [0, 1, 5], [0, 5, 4],
        [3, 1, 2], [3, 5, 1], [3, 2, 4], [3, 4, 5],
    ],
)

# Cube with vertices at (+-1, +-1, +-1).
cube = Polyhedron(
    vertex_positions=[
        (-1, -1, -1), (-1, -1, +1), (-1, +1, -1), (-1, +1, +1),
        (+1, -1, -1), (+1, -1, +1), (+1, +1, -1), (+1, +1, +1),
    ],
    triangles=[
        [1, 3, 2], [1, 0, 4], [1, 5, 7],
        [2, 0, 1], [2, 6, 4], [2, 3, 7],
        [4, 5, 1], [4, 0, 2], [4, 6, 7],
        [7, 3, 1], [7, 6, 2], [7, 5, 4],
    ],
)

# Pair of cubes, sharing a common vertex at the origin.
pair_of_cubes = Polyhedron(
    vertex_positions=[
        (-1, -1, -1), (-1, -1, 0), (-1, 0, -1), (-1, 0, 0),
        (0, -1, -1), (0, -1, 0), (0, 0, -1), (0, 0, 0),
        (0, 0, 1), (0, 1, 0), (0, 1, 1),
        (1, 0, 0), (1, 0, 1), (1, 1, 0), (1, 1, 1),
    ],
    triangles=[
        [1, 3, 2], [1, 0, 4], [1, 5, 7],
        [2, 0, 1], [2, 6, 4], [2, 3, 7],
        [4, 5, 1], [4, 0, 2], [4, 6, 7],
        [7, 3, 1], [7, 6, 2], [7, 5, 4],
        [8, 10, 9], [8, 7, 11], [8, 12, 14],
        [9, 7, 8], [9, 13, 11], [9, 10, 14],
        [11, 12, 8], [11, 7, 9], [11, 13, 14],
        [14, 10, 8], [14, 13, 9], [14, 12, 11],
    ],
)

# Two stacked cuboids, one directly above the other.

aligned_stacked_cuboids = Polyhedron(
    vertex_positions=[
        (0, 0, 0), (0, 0, 1), (0, 3, 0), (0, 3, 1),
        (3, 0, 0), (3, 0, 1), (3, 3, 0), (3, 3, 1),
        (0, 0, 2), (0, 0, 3), (0, 3, 2), (0, 3, 3),
        (3, 0, 2), (3, 0, 3), (3, 3, 2), (3, 3, 3),
    ],
    triangles=(
        cube.triangles +
        [[x+8, y+8, z+8] for x, y, z in cube.triangles]
    ),
)

# Similar, but with the cuboids misaligned.

misaligned_stacked_cuboids = Polyhedron(
    vertex_positions=[
        (0, 0, 0), (0, 0, 1), (0, 2, 0), (0, 2, 1),
        (2, 0, 0), (2, 0, 1), (2, 2, 0), (2, 2, 1),
        (1, 1, 2), (1, 1, 3), (1, 3, 2), (1, 3, 3),
        (3, 1, 2), (3, 1, 3), (3, 3, 2), (3, 3, 3),
    ],
    triangles=(
        cube.triangles +
        [[x+8, y+8, z+8] for x, y, z in cube.triangles]
    ),
)

# Hollow cube: surface consists of two cubes, one facing outwards
# and one facing inwards.
hollow_cube = Polyhedron(
    vertex_positions=[
        (0, 0, 0), (0, 0, 3), (0, 3, 0), (0, 3, 3),
        (3, 0, 0), (3, 0, 3), (3, 3, 0), (3, 3, 3),
        (1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2),
        (2, 1, 1), (2, 1, 2), (2, 2, 1), (2, 2, 2),
    ],
    triangles=(
        # Outer cube.
        cube.triangles +
        # Inner cube: reverse orientation.
        [[z+8, y+8, x+8] for x, y, z in cube.triangles]
    ),
)

# Torus.
torus = Polyhedron(
    vertex_positions=[
        # Outer square, bottom face.
        (0, 0, 0), (0, 3, 0), (3, 0, 0), (3, 3, 0),
        # Inner square, bottom face.
        (1, 1, 0), (1, 2, 0), (2, 1, 0), (2, 2, 0),
        # Outer square, top face.
        (0, 0, 1), (0, 3, 1), (3, 0, 1), (3, 3, 1),
        # Inner square, top face.
        (1, 1, 1), (1, 2, 1), (2, 1, 1), (2, 2, 1),
    ],
    triangles=[
        # Top face.
        [8, 14, 12], [14, 8, 10], [10, 15, 14], [15, 10, 11],
        [11, 13, 15], [13, 11, 9], [9, 12, 13], [12, 9, 8],
        # Bottom face.
        [4, 6, 0], [2, 0, 6], [6, 7, 2], [3, 2, 7],
        [7, 5, 3], [1, 3, 5], [5, 4, 1], [0, 1, 4],
        # Outer faces.
        [0, 2, 10], [10, 8, 0], [2, 3, 11], [11, 10, 2],
        [3, 1, 9], [9, 11, 3], [1, 0, 8], [8, 9, 1],
        # Inner faces.
        [4, 12, 14], [14, 6, 4], [6, 14, 15], [15, 7, 6],
        [7, 15, 13], [13, 5, 7], [5, 13, 12], [12, 4, 5],
    ],
)


class TestPolyhedron(unittest.TestCase):
    def test_tetrahedron(self):
        xs = ys = zs = [0.25 * v for v in range(-1, 6)]
        points = [(x, y, z) for x in xs for y in ys for z in zs]
        for point in points:
            class_ = tetrahedron_classify(point)
            if class_ == "inside":
                self.assertEqual(tetrahedron.winding_number(point), 1)
            elif class_ == "outside":
                self.assertEqual(tetrahedron.winding_number(point), 0)
            elif class_ == "boundary":
                # Point is on the boundary.
                with self.assertRaises(ValueError):
                    tetrahedron.winding_number(point)
            else:
                assert False, "should never get here"

    def test_cube(self):
        # Check volume
        self.assertEqual(cube.volume(), 8)

        def classify(point):
            x, y, z = point
            if -1 < x < 1 and -1 < y < 1 and -1 < z < 1:
                return "inside"
            if -1 <= x <= 1 and -1 <= y <= 1 and -1 <= z <= 1:
                return "boundary"
            return "outside"

        # Quarter-integer boundaries from -1.25 to 1.25 inclusive.
        xs = ys = zs = [0.25 * v for v in range(-5, 6)]
        points = [(x, y, z) for x in xs for y in ys for z in zs]
        for point in points:
            class_ = classify(point)
            if class_ == "inside":
                self.assertEqual(cube.winding_number(point), 1)
            elif class_ == "outside":
                self.assertEqual(cube.winding_number(point), 0)
            elif class_ == "boundary":
                # Point is on the boundary.
                with self.assertRaises(ValueError):
                    cube.winding_number(point)
            else:
                assert False, "should never get here"

    def test_octahedron(self):
        self.assertEqual(octahedron.volume(), 4.0 / 3.0)

        def classify(point):
            x, y, z = point
            s = abs(x) + abs(y) + abs(z)
            if s < 1:
                return "inside"
            elif s == 1:
                return "boundary"
            else:
                return "outside"

        # Quarter-integer boundaries from -1.25 to 1.25 inclusive.
        xs = ys = zs = [0.25 * v for v in range(-5, 6)]
        points = [(x, y, z) for x in xs for y in ys for z in zs]
        for point in points:
            class_ = classify(point)
            if class_ == "inside":
                self.assertEqual(octahedron.winding_number(point), 1)
            elif class_ == "outside":
                self.assertEqual(octahedron.winding_number(point), 0)
            elif class_ == "boundary":
                # Point is on the boundary.
                with self.assertRaises(ValueError):
                    octahedron.winding_number(point)
            else:
                assert False, "never get here"

    def test_pair_of_cubes(self):
        self.assertEqual(pair_of_cubes.volume(), 2.0)

        def classify(point):
            x, y, z = point
            if -1 < x < 0 and -1 < y < 0 and -1 < z < 0:
                return "inside"
            if 0 < x < 1 and 0 < y < 1 and 0 < z < 1:
                return "inside"
            if -1 <= x <= 0 and -1 <= y <= 0 and -1 <= z <= 0:
                return "boundary"
            if 0 <= x <= 1 and 0 <= y <= 1 and 0 <= z <= 1:
                return "boundary"
            return "outside"

        # Quarter-integer boundaries from -1.25 to 1.25 inclusive.
        xs = ys = zs = [0.25 * v for v in range(-5, 6)]
        points = [(x, y, z) for x in xs for y in ys for z in zs]
        for point in points:
            class_ = classify(point)
            if class_ == "inside":
                self.assertEqual(pair_of_cubes.winding_number(point), 1)
            elif class_ == "boundary":
                with self.assertRaises(ValueError):
                    pair_of_cubes.winding_number(point)
            elif class_ == "outside":
                self.assertEqual(pair_of_cubes.winding_number(point), 0)
            else:
                assert False, "never get here"

    def test_hollow_cube(self):
        self.assertEqual(hollow_cube.volume(), 26.0)

        def classify(point):
            x, y, z = point
            if 1 < x < 2 and 1 < y < 2 and 1 < z < 2:
                return "outside"
            if 1 <= x <= 2 and 1 <= y <= 2 and 1 <= z <= 2:
                return "boundary"
            if 0 < x < 3 and 0 < y < 3 and 0 < z < 3:
                return "inside"
            if 0 <= x <= 3 and 0 <= y <= 3 and 0 <= z <= 3:
                return "boundary"
            return "outside"

        xs = ys = zs = [0.25 * v for v in range(-1, 14)]
        points = [(x, y, z) for x in xs for y in ys for z in zs]
        for point in points:
            class_ = classify(point)
            if class_ == "inside":
                self.assertEqual(hollow_cube.winding_number(point), 1)
            elif class_ == "boundary":
                with self.assertRaises(ValueError):
                    hollow_cube.winding_number(point)
            elif class_ == "outside":
                self.assertEqual(hollow_cube.winding_number(point), 0)
            else:
                assert False, "never get here"

    def test_torus(self):
        self.assertEqual(torus.volume(), 8.0)

        def classify(point):
            x, y, z = point
            if 0 < z < 1 and (0 < x < 1 or 2 < x < 3) and 0 < y < 3:
                return "inside"
            if 0 < z < 1 and (0 < y < 1 or 2 < y < 3) and 0 < x < 3:
                return "inside"
            if 0 <= z <= 1 and (0 <= x <= 1 or 2 <= x <= 3) and 0 <= y <= 3:
                return "boundary"
            if 0 <= z <= 1 and (0 <= y <= 1 or 2 <= y <= 3) and 0 <= x <= 3:
                return "boundary"
            return "outside"

        xs = ys = [0.25 * v for v in range(-1, 14)]
        zs = [0.25 * v for v in range(-1, 6)]
        points = [(x, y, z) for x in xs for y in ys for z in zs]
        for point in points:
            class_ = classify(point)
            if class_ == "inside":
                self.assertEqual(torus.winding_number(point), 1)
            elif class_ == "boundary":
                with self.assertRaises(ValueError):
                    torus.winding_number(point)
            elif class_ == "outside":
                self.assertEqual(torus.winding_number(point), 0)
            else:
                assert False, "never get here"

    def test_aligned_stacked_cuboids(self):
        self.assertEqual(aligned_stacked_cuboids.volume(), 18.0)

        def classify(point):
            x, y, z = point
            if 0 < x < 3 and 0 < y < 3 and 0 < z < 1:
                return "inside"
            if 0 < x < 3 and 0 < y < 3 and 2 < z < 3:
                return "inside"
            if 0 <= x <= 3 and 0 <= y <= 3 and 0 <= z <= 1:
                return "boundary"
            if 0 <= x <= 3 and 0 <= y <= 3 and 2 <= z <= 3:
                return "boundary"
            return "outside"

        xs = ys = zs = [0.25 * v for v in range(-1, 14)]
        points = [(x, y, z) for x in xs for y in ys for z in zs]
        for point in points:
            class_ = classify(point)
            if class_ == "inside":
                self.assertEqual(
                    aligned_stacked_cuboids.winding_number(point), 1)
            elif class_ == "boundary":
                with self.assertRaises(ValueError):
                    aligned_stacked_cuboids.winding_number(point)
            elif class_ == "outside":
                self.assertEqual(
                    aligned_stacked_cuboids.winding_number(point), 0)
            else:
                assert False, "never get here"

    def test_misaligned_stacked_cuboids(self):
        self.assertEqual(misaligned_stacked_cuboids.volume(), 8.0)

        def classify(point):
            x, y, z = point
            if 0 < x < 2 and 0 < y < 2 and 0 < z < 1:
                return "inside"
            if 1 < x < 3 and 1 < y < 3 and 2 < z < 3:
                return "inside"
            if 0 <= x <= 2 and 0 <= y <= 2 and 0 <= z <= 1:
                return "boundary"
            if 1 <= x <= 3 and 1 <= y <= 3 and 2 <= z <= 3:
                return "boundary"
            return "outside"

        xs = ys = zs = [0.25 * v for v in range(-1, 14)]
        points = [(x, y, z) for x in xs for y in ys for z in zs]
        for point in points:
            class_ = classify(point)
            if class_ == "inside":
                self.assertEqual(
                    misaligned_stacked_cuboids.winding_number(point), 1)
            elif class_ == "boundary":
                with self.assertRaises(ValueError):
                    misaligned_stacked_cuboids.winding_number(point)
            elif class_ == "outside":
                self.assertEqual(
                    misaligned_stacked_cuboids.winding_number(point), 0)
            else:
                assert False, "never get here"


if __name__ == '__main__':
    unittest.main()
