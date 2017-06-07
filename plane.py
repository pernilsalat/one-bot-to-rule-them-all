import math
import numpy

def inted(matrix):
    return list(map(list,map(lambda x: map(int,x),matrix)))

class Plane:
    def __init__(self):
        self.origin = (0,0)
        self.Transform = numpy.identity(3)

    def moveTo(self, point):
        p  = self.getCoord(point)
        x, y = p[0][0], p[1][0]
        Tran = numpy.identity(3)
        Tran[0][2] = -x
        Tran[1][2] = -y
        self.Transform = Tran @ self.Transform

    def lookAt(self, point):
        p = self.getCoord(point)
        x, y = p[0][0], p[1][0]
        slope = -float(x) / float(y)
        angle = math.atan(slope)
        Tran = numpy.zeros((3,3))
        Tran[0][0] = math.cos(angle)
        Tran[0][1] = math.sin(angle)
        Tran[1][0] = -Tran[0][1]
        Tran[1][1] = Tran[0][0]
        Tran[2][2] = 1
        self.Transform = Tran @ self.Transform

    def lookAt18(self, point):
        p = self.getCoord(point)
        x, y = p[0][0], p[1][0]
        slope = -float(x) / float(y)
        angle = math.atan(slope)
        if angle > 18:
            angle = 18
        elif angle < -18:
            angle = -18
        Tran = numpy.zeros((3,3))
        Tran[0][0] = math.cos(angle)
        Tran[0][1] = math.sin(angle)
        Tran[1][0] = -Tran[0][1]
        Tran[1][1] = Tran[0][0]
        Tran[2][2] = 1
        self.Transform = Tran @ self.Transform

    def getCoord(self, point):
        x, y = point
        return self.Transform @ numpy.matrix([[x],[y],[1]])

    def realCoord(self, point):
        x, y = point
        inv = numpy.linalg.inv(self.Transform) 
        return inv @ numpy.matrix([[x],[y],[1]])
