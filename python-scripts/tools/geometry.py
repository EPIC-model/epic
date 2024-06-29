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
