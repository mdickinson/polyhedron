"""
Class representing a Polyhedron, described by a triangular mesh.

"""

def sign(x):
    """
    Return 1 if x is positive, -1 if it's negative, and 0 if it's zero.

    """
    return (x > 0) - (x < 0)


def triangle_chain(v1, v2, v3, origin):
    """
    Map the triangle v1-v2-v3 to the corresponding 2-chain in cellular
    homology.

    Return coefficient of [U] in the result.

    Raise ValueError if the face contains the origin.

    """
    # Result depends on various minors of the matrix:
    #
    # ( v1[0] v2[0] v3[0] origin[0] )
    # ( v1[1] v2[1] v3[1] origin[1] )
    # ( v1[2] v2[2] v3[2] origin[2] )
    # (   1     1     1      1      )

    m1_0 = v1[0] - origin[0]
    m1_1 = v1[1] - origin[1]
    m1_2 = v1[2] - origin[2]
    if m1_0 == m1_1 == m1_2 == 0:
        raise ValueError("v1 == origin")
    v1positive = m1_0 > 0 or m1_0 == 0 and (m1_1 > 0 or m1_1 == 0 and m1_2 > 0)

    m2_0 = v2[0] - origin[0]
    m2_1 = v2[1] - origin[1]
    m2_2 = v2[2] - origin[2]
    if m2_0 == m2_1 == m2_2 == 0:
        raise ValueError("v2 == origin")
    v2positive = m2_0 > 0 or m2_0 == 0 and (m2_1 > 0 or m2_1 == 0 and m2_2 > 0)

    m3_0 = v3[0] - origin[0]
    m3_1 = v3[1] - origin[1]
    m3_2 = v3[2] - origin[2]
    if m3_0 == m3_1 == m3_2 == 0:
        raise ValueError("v3 == origin")
    v3positive = m3_0 > 0 or m3_0 == 0 and (m3_1 > 0 or m3_1 == 0 and m3_2 > 0)

    if v1positive == v2positive == v3positive:
        return 0

    # Look for positive edges from N to P,
    # and negative edges from P to N.
    face_boundary = 0

    edge_boundary = v2positive - v1positive
    m12_01 = d1 = m1_0 * m2_1 - m1_1 * m2_0
    if edge_boundary:
        # More minors: m12_01, m12_02, m12_12
        m12_02 = d2 = m1_0 * m2_2 - m1_2 * m2_0
        m12_12 = d3 = m1_1 * m2_2 - m1_2 * m2_1
        if m12_01 == m12_02 == m12_12 == 0:
            raise ValueError("v1, v2, origin collinear")
        e12_positive = d1 < 0 or d1 == 0 and (d2 < 0 or d2 == 0 and d3 < 0)
        if e12_positive == v2positive:
            face_boundary += edge_boundary

    edge_boundary = v3positive - v2positive
    m23_01 = d1 = m2_0 * m3_1 - m2_1 * m3_0
    if edge_boundary:
        m23_02 = d2 = m2_0 * m3_2 - m2_2 * m3_0
        m23_12 = d3 = m2_1 * m3_2 - m2_2 * m3_1
        if d1 == d2 == d3 == 0:
            raise ValueError("v2, v3, origin collinear")
        e23_positive = d1 < 0 or d1 == 0 and (d2 < 0 or d2 == 0 and d3 < 0)
        if e23_positive == v3positive:
            face_boundary += edge_boundary

    edge_boundary = v1positive - v3positive
    m31_01 = d1 = m3_0 * m1_1 - m3_1 * m1_0
    if edge_boundary:
        m31_02 = d2 = m3_0 * m1_2 - m3_2 * m1_0
        m31_12 = d3 = m3_1 * m1_2 - m3_2 * m1_1
        if d1 == d2 == d3 == 0:
            raise ValueError("v3, v1, origin collinear")
        e31_positive = d1 < 0 or d1 == 0 and (d2 < 0 or d2 == 0 and d3 < 0)
        if e31_positive == v1positive:
            face_boundary += edge_boundary

    if not face_boundary:
        return 0

    ts = m12_01 * m3_2 + m23_01 * m1_2 + m31_01 * m2_2
    if not ts:
        raise ValueError("v1, v2, v3, origin are coplanar")
    if (face_boundary > 0) == (ts > 0):
        return face_boundary
    else:
        return 0


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
            triangle_chain(v1, v2, v3, point)
            for v1, v2, v3 in self.triangle_positions())
