"""Robust point-in-polyhedron test.

Given an closed, oriented surface in R^3 described by a triangular mesh, the
code below gives a robust algorithm for determining whether a given point is
inside, on the boundary of, or outside, the surface.  The algorithm should give
correct results even in degenerate cases, and applies to disconnected
polyhedra, non simply-connected surfaces, and so on.  There are no requirements
for the surface to be convex, simple, connected or simply-connected.

More precisely, we give a method for computing the *winding number* of a closed
oriented surface around a point P.  Roughly speaking, the winding number of the
closed oriented surface S around a point P not on S is the number of times that
the surface encloses that point; for a simple outward-oriented surface (like
that of a convex polyhedron, for example), the winding number will be 1 for
points inside the surface and 0 for points outside.

For a precise definition of winding number, we can turn to algebraic topology:
our oriented surface is presented as a collection of combinatorial data
defining abstract vertices, edges and triangles, together with a mapping of
those vertices to R^3.  The combinatorial data describe a simplicial complex C,
and assuming that P doesn't lie on the surface, the mapping of the vertices to
R^3 gives a continuous map from the geometric realization of C to R^3 - {P}.
This in turn induces a map on second homology groups:

   H^2(C, Z) -> H^2(R^3 - {P}, Z)

and by taking the usual right-handed orientation in R^3 we identify H^2(R^3 -
{P}, Z) with Z.  The image of [S] under this map gives the winding number.  In
particular, the well-definedness of the winding number does not depend on
topological properties of the embedding: it doesn't matter if the surface is
self-intersecting, or has degenerate triangles.  The only condition is that P
does not lie on the surface S.

Algorithm
---------
The algorithm is based around the usual method of ray-casting: we take a
vertical line L through P and count the intersections of this line with the
triangles of the surface, keeping track of orientations as we go.  Let's ignore
corner cases for a moment and assume that:

(1) P does not lie on the surface, and
(2) for each triangle T (thought of as a closed subset of R^3) touched by
    our vertical line L, L meets the interior of T in exactly one point Q

Then there are four possibilities for each such triangle T:

1. T lies *above* P and is oriented *upwards* (*away* from P).
2. T lies *above* P and is oriented *downwards* (*towards* P).
3. T lies *below* P and is oriented *downwards* (*away* from P).
4. T lies *below* P and is oriented *upwards* (*towards* P).

Let's write N1, N2, N3 and N4 for the counts of triangles satisfying conditions
1, 2, 3 and 4 respectively.  Since we have a closed surface, these numbers
are not independent; they satisfy the relation:

    N1 + N4 == N2 + N3

That is, the number of upward-facing triangles must match the number of
downward-facing triangles.  The winding number w is then given by:

    w = N1 - N2 == N3 - N4

In the code below, we simply compute 2*w = (N1 + N3) - (N2 + N4), so each
triangle oriented away from P contributes 1 to 2w, while each triangle oriented
towards P contributes -1.

Determining whether a triangle QRS points away from or towards P involves
finding the sign of the determinant of the 4x4 matrix

   | Qx Rx Sx Px |
   | Qy Ry Sy Py |
   | Qz Rz Sz Pz |
   |  1  1  1  1 |

or equivalently of the 3x3 matrix

   | Qx-Px Rx-Px Sx-Px |
   | Qy-Py Ry-Py Sy-Py |
   | Qz-Pz Rz-Pz Sz-Pz |


Making the algorithm robust
---------------------------
Now we describe how to augment the basic algorithm described above to include:

- correct treatment of corner cases (vertical triangles, cases where L meets an
  edge or vertex directly, etc.)

- detection of the case where the point lies directly on the surface.

It turns out that to make the algorithm robust, all we need to do is be careful
and consistent about classifying vertices, edges and triangles.  We do this as
follows:

- Each vertex of the surface that's not equal to P is considered *positive* if
  its coordinates are lexicographically greater than P, and *negative*
  otherwise.

- Each edge QR of the surface that's not collinear with P is considered
  *positive* if its projection to the x-y plane passes to the right of P (i.e.,
  if Q-R-P represents a counterclockwise turn), and *negative* if it passes to
  the left.  If the projected edge passes *through* P then the actual edge must
  intersect L; in that case we consider the edge to be *positive* if it
  intersects L at a point *below* P, and *negative* if it intersects L at a
  point *above* P.

- Each triangle QRS of the surface that's not coplanar with P is considered
  *positive* if its normal points away from P, and *negative* if its normal
  points towards P.

Now to compute the contribution of any given triangle to the total winding
number:

1. Classify the vertices of the triangle.  At the same time, we can check that
   none of the vertices is equal to P.  If all vertices have the same sign,
   then the winding number contribution is zero.

2. Assuming that the vertices do not all have the same sign, two of the three
   edges connect two differently-signed vertices.  Classify both those edges
   (and simultaneously check that they don't pass through P).  If the edges
   have opposite classification, then the winding number contribution is zero.

3. Now two of the edges have the same sign: classify the triangle itself.  If
   the triangle is positive it contributes 1/2 to the winding number total; if
   negative it contributes -1/2.  In practice we count contributions of 1 and
   -1, and halve the total at the end.

Note that an edge between two like-signed vertices can never pass through P, so
there's no need to check the third edge in step 2.  Similarly, a triangle whose
edge-cycle is trivial can't contain P in its interior.

"""


def sign(x):
    """
    Return 1 if x is positive, -1 if it's negative, and 0 if it's zero.

    """
    return (x > 0) - (x < 0)


def vertex_sign(v1, origin):
    result = (
        sign(v1[0] - origin[0]) or
        sign(v1[1] - origin[1]) or
        sign(v1[2] - origin[2]))
    if not result:
        raise ValueError("vertex coincides with origin")
    return result


def edge_sign(v1, v2, origin):
    result = (
        sign((v1[0] - origin[0]) * (v2[1] - origin[1]) -
             (v1[1] - origin[1]) * (v2[0] - origin[0])) or
        sign((v1[0] - origin[0]) * (v2[2] - origin[2]) -
             (v1[2] - origin[2]) * (v2[0] - origin[0])) or
        sign((v1[1] - origin[1]) * (v2[2] - origin[2]) -
             (v1[2] - origin[2]) * (v2[1] - origin[1])))
    if not result:
        raise ValueError("vertices collinear with origin")
    return result


def triangle_sign(v1, v2, v3, origin):
    m1_0 = v1[0] - origin[0]
    m1_1 = v1[1] - origin[1]
    m1_2 = v1[2] - origin[2]
    m2_0 = v2[0] - origin[0]
    m2_1 = v2[1] - origin[1]
    m2_2 = v2[2] - origin[2]
    m3_0 = v3[0] - origin[0]
    m3_1 = v3[1] - origin[1]
    m3_2 = v3[2] - origin[2]
    m12_01 = m1_0 * m2_1 - m1_1 * m2_0
    m23_01 = m2_0 * m3_1 - m2_1 * m3_0
    m31_01 = m3_0 * m1_1 - m3_1 * m1_0
    result = sign(m12_01 * m3_2 + m23_01 * m1_2 + m31_01 * m2_2)
    if not result:
        raise ValueError("vertices coplanar with origin")
    return result


def triangle_chain(v1, v2, v3, origin):
    """
    Return the contribution of this triangle to the winding number.

    Raise ValueError if the face contains the origin.

    """
    v1sign = vertex_sign(v1, origin)
    v2sign = vertex_sign(v2, origin)
    v3sign = vertex_sign(v3, origin)

    face_boundary = 0
    if v1sign != v2sign:
        face_boundary += edge_sign(v1, v2, origin)
    if v2sign != v3sign:
        face_boundary += edge_sign(v2, v3, origin)
    if v3sign != v1sign:
        face_boundary += edge_sign(v3, v1, origin)
    if not face_boundary:
        return 0

    return triangle_sign(v1, v2, v3, origin)


class Polyhedron(object):
    def __init__(self, triangles, vertex_positions):
        # Vertex positions in R^3.
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

    def volume(self):
        """
        Return the volume of this polyhedron.

        """
        acc = 0
        for p1, p2, p3 in self.triangle_positions():
            # Twice the area of the projection onto the x-y plane.
            det = ((p2[1] - p3[1]) * (p1[0] - p3[0]) -
                   (p2[0] - p3[0]) * (p1[1] - p3[1]))
            # Three times the average height.
            height = p1[2] + p2[2] + p3[2]
            acc += det * height
        return acc / 6.0

    def winding_number(self, point):
        """Determine the winding number of *self* around the given point.

        """
        return sum(
            triangle_chain(v1, v2, v3, point)
            for v1, v2, v3 in self.triangle_positions()) // 2
