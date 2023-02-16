import numpy as np

class plane:
    def __init__(self, normal, point):
        self.normal = np.asarray(normal)
        self.normal = self.normal/ np.linalg.norm(self.normal, 2)
        self.point = np.asarray(point)
        d = - self.point.dot(self.normal)
        self.plane = np.array([self.normal[0], self.normal[1], self.normal[2], d])
        self.name = 'arbitrary-plane'
        #print("Plane: ax + by + cz + d = 0")
        #print('a:', self.plane[0])
        #print('b:', self.plane[1])
        #print('c:', self.plane[2])
        #print('d:', self.plane[3])


class xy_plane(plane):
    def __init__(self, z=0.0):
        super(xy_plane, self).__init__(normal=[0, 0, 1], point=[0, 0, z])
        self.name = 'xy-plane'
        self.loc = z


class xz_plane(plane):
    def __init__(self, y=0.0):
        super(xz_plane, self).__init__(normal=[0, 1, 0], point=[0, y, 0])
        self.name = 'xz-plane'
        self.loc = y


class yz_plane(plane):
    def __init__(self, x=0.0):
        super(yz_plane, self).__init__(normal=[1, 0, 0], point=[x, 0, 0])
        self.name = 'yz-plane'
        self.loc = x

class ellipsoid:
    def __init__(self, xc, A):
        A = np.asarray(A)

        if not A.shape == (3, 3):
            raise ValueError('Not a 3x3 matrix')
        #if not (A.transpose() == A).all():
            #raise ValueError('Not a symmetriy matrix')

        # implicit: Ax**2 + Bxy + Cxz + Dy**2 + Eyz + Fz**2 + Gx + Hy + Iz + J = 0
        self.A = A[0, 0]
        self.B = 2.0 * A[0, 1]
        self.C = 2.0 * A[0, 2]
        self.D = A[1, 1]
        self.E = 2.0 * A[1, 2]
        self.F = A[2, 2]
        self.G = -2.0 * (A[0, 0] * xc[0] + A[0, 1] * xc[1] + A[0, 2] * xc[2])
        self.H = -2.0 * (A[0, 1] * xc[0] + A[1, 1] * xc[1] + A[1, 2] * xc[2])
        self.I = -2.0 * (A[0, 2] * xc[0] + A[1, 2] * xc[1] + A[2, 2] * xc[2])
        self.J = 2.0 * A[0, 1] * xc[0] * xc[1] + A[0, 0] * xc[0]**2 \
               + 2.0 * A[0, 2] * xc[0] * xc[2] + A[1, 1] * xc[1]**2 \
               + 2.0 * A[1, 2] * xc[1] * xc[2] + A[2, 2] * xc[2]**2 \
               - 1.0
        #print("Ellipsoid: Ax**2 + Bxy + Cxz + Dy**2 + Eyz + Fz**2 + Gx + Hy + Iz + J = 0")
        #print('A:', self.A)
        #print('B:', self.B)
        #print('C:', self.C)
        #print('D:', self.D)
        #print('E:', self.E)
        #print('F:', self.F)
        #print('G:', self.G)
        #print('H:', self.H)
        #print('I:', self.I)
        #print('J:', self.J)

    def intersect(self, pl):
        a = pl.plane[0]
        b = pl.plane[1]
        c = pl.plane[2]
        d = pl.plane[3]
        A = self.A
        B = self.B
        C = self.C
        D = self.D
        E = self.E
        F = self.F
        G = self.G
        H = self.H
        I = self.I
        J = self.J

        ell = ellipse()
        if not c == 0 and pl.name == 'xy-plane':
            ell.A = A - C * a / c + F * a**2 / c**2
            ell.B = B - C * b / c - E * a / c + 2.0 * a * b * F / c**2
            ell.C = D - E * b / c + F * b**2 / c**2
            ell.D = - C * d / c + 2.0 * d * a * F / c**2 + G - I * a / c
            ell.E = - E * d / c + 2.0 * d * b * F / c**2 + H - I * b / c
            ell.F = - I * d / c + J + F * d**2 / c**2
        elif not b == 0 and pl.name == 'xz-plane':
            ell.A = A - B * a / b + D * a**2 / b**2
            ell.B = C - B * c / b - E * a / b + 2.0 * a * c * D / b**2
            ell.C = F - E * c / b + D * c**2 / b**2
            ell.D = - B * d / b + 2.0 * d * a * D / b**2 + G - H * a / b
            ell.E = - E * d / b + 2.0 * d * c * D / b**2 + I - H * c / b
            ell.F = - H * d / b + J + D * d**2 / b**2
        elif not a == 0 and pl.name == 'yz-plane':
            ell.A = D - B * b / a + A * b**2 / a**2
            ell.B = E - B * c / a - C * b / a + 2.0 * b * c * A / a**2
            ell.C = F - C * c / a + A * c**2 / a**2
            ell.D = - B * d / a + 2.0 * b * d * A / a**2 + H - G * b / a
            ell.E = - C * d / a + 2.0 * c * d * A / a**2 + I - G * c / a
            ell.F = - G * d / a + J + A * d**2 / a**2
        else:
            print("Not yet implemented.")
            exit()
        ell.axis = pl.name
        return ell

# 2 Jan. 2023
# https://en.wikipedia.org/wiki/Ellipse#General_ellipse
class ellipse:

    def __init__(self):
        # Ax**2 + Bxy + Cy**2 + Dx + Ey + F = 0
        self.A = 0.0
        self.B = 0.0
        self.C = 0.0
        self.D = 0.0
        self.E = 0.0
        self.F = 0.0
        self.axis = None

    def is_real(self):
        det = self._determinant()
        return self.C * det < 0.0

    def _determinant(self):
        return (self.A * self.C - 0.25 * self.B**2) * self.F \
            + 0.25 * self.B * self.E * self.D \
            - 0.25 * self.C * self.D * self.D \
            - 0.25 * self.A * self.E * self.E

    @property
    def centre(self):
        denom = self.B**2 - 4.0 * self.A * self.C
        xo = (2.0 * self.C * self.D - self.B * self.E) / denom
        yo = (2.0 * self.A * self.E - self.B * self.D) / denom
        return (xo, yo)

    @property
    def semi_axes(self):
        denom = self.B**2 - 4.0 * self.A * self.C
        A = self.A * self.E**2 + self.C * self.D**2 - self.B * self.D * self.E + denom * self.F
        B = self.A + self.C + np.sqrt((self.A - self.C)**2 + self.B**2)
        a = - np.sqrt(2.0 * A * B) / denom
        B = self.A + self.C - np.sqrt((self.A - self.C)**2 + self.B**2)
        b = - np.sqrt(2.0 * A * B) / denom
        return (a, b)

    @property
    def angle(self):
        A = self.A
        B = self.B
        C = self.C
        if B == 0.0 and A < C:
            phi = 0.0
        elif B == 0.0 and A > C:
            phi = np.pi/2.0
        elif not B == 0.0:
            phi = np.arctan2(C - A - np.sqrt((A - C)**2 + B**2), B)
        else:
            raise ValueError('Cannot evaluate ellipse angle.')

        if self.axis == 'xz-plane':
            phi = np.pi - phi
        elif self.axis == 'yz-plane':
            phi = np.pi/2.0 + phi
        return phi
