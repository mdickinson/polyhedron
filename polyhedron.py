"""
Class representing a Polyhedron, described by a triangular mesh.

"""


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

    Return coefficient of B in C_0 = Z + Z.

    Raise ValueError if v1 == origin.

    """
    d = (sign(v1[0] - origin[0]) or
         sign(v1[1] - origin[1]) or sign(v1[2] - origin[2]))
    if not d:
        raise ValueError("Vertex contains origin")
    return 0 if d > 0 else 1


def edge_chain(v1, v2, origin):
    """
    Map the edge v1-v2 to the corresponding 1-chain in cellular homology.

    Return coefficient of L in C_1 = Z + Z.

    Raise ValueError if the edge goes through the origin.

    """
    edge_boundary = vertex_chain(v2, origin) - vertex_chain(v1, origin)
    if not edge_boundary:
        return 0

    d = (ccw2(v1, v2, origin) or
         ccw2((v1[0], v1[2]), (v2[0], v2[2]), (origin[0], origin[2])) or
         ccw2((v1[1], v1[2]), (v2[1], v2[2]), (origin[1], origin[2])))
    if not d:
        raise ValueError("Edge contains origin")

    return 0 if d * edge_boundary > 0 else d


def face_chain(v1, v2, v3, origin):
    """
    Map the triangle v1-v2-v3 to the corresponding 2-chain in cellular
    homology.

    Return coefficient of D in the result.

    Raise ValueError if the face contains the origin.

    """
    face_boundary = sum(
        edge_chain(start, end, origin)
        for start, end in [(v1, v2), (v2, v3), (v3, v1)])

    if not face_boundary:
        return 0

    d = ccw3(v1, v2, v3, origin)
    if not d:
        raise ValueError("Face contains origin")

    return 0 if d * face_boundary > 0 else d


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
        """Determine whether the given point is inside the polyhedron.

        More precisely, determine the winding number of the polyhedron around
        the point.

        """
        return sum(
            face_chain(v1, v2, v3, point)
            for v1, v2, v3 in self.triangle_positions())
