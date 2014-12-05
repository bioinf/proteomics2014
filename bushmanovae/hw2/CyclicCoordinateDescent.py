__author__ = 'letovesnoi'

import sys

def main():
    inputPDB = sys.argv[1]
    outputPDB = sys.argv[2]
    coordinates, types = getXYZTypeFromPDBline(inputPDB)
    lastCoordinate = sys.argv[3].split(' ')
    last_point = Point(float(lastCoordinate[0]), float(lastCoordinate[1]), float(lastCoordinate[2]))
    newCoordinates = getNewConfiguration(coordinates, types, last_point)
    printPDBResults(inputPDB, outputPDB, newCoordinates)

def getXYZTypeFromPDBline(inputPDB):
    coordinates = []
    types = []
    with open(inputPDB, 'r') as fin:
        for line in fin:
            if line[0:6] == 'ATOM  ':
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                type = line[77:80].strip()
                coordinates.append((x, y, z))
                types.append(type)
    return coordinates, types

def printPDBResults(inputPDB, outputPDB, newCoordinates):
    iAtom = 0
    with open(inputPDB, 'r') as fin:
        with open(outputPDB, 'w') as fout:
            for line in fin:
                if line[0:6] == 'ATOM  ':
                    fout.write(line[:30] + '{0:>8.3f}'.format(newCoordinates[iAtom][0]) +
                               '{0:>8.3f}'.format(newCoordinates[iAtom][1]) +
                               '{0:>8.3f}'.format(newCoordinates[iAtom][2]) + line[54:])
                    iAtom += 1
                else:
                    fout.write(line)

class Point(object):
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def rotate(self, point0, direction_vector, teta):
        import math
        x = (point0.x * (direction_vector.y * direction_vector.y + direction_vector.z * direction_vector.z) -
             direction_vector.x * (point0.y * direction_vector.y + point0.z * direction_vector.z - direction_vector.x * self.x
                                   - direction_vector.y * self.y - direction_vector.z * self.z)) * (1 - math.cos(teta)) + self.x * math.cos(teta) + \
            (-point0.z * direction_vector.y + point0.y * direction_vector.z - direction_vector.z * self.y + direction_vector.y * self.z) * math.sin(teta)

        y = (point0.y * (direction_vector.x * direction_vector.x + direction_vector.z * direction_vector.z) -
             direction_vector.y * (point0.x * direction_vector.x + point0.z * direction_vector.z - direction_vector.x * self.x
                                   - direction_vector.y * self.y - direction_vector.z * self.z)) * (1 - math.cos(teta)) + self.y * math.cos(teta) + \
            (point0.z * direction_vector.x + point0.x * direction_vector.z - direction_vector.z * self.x - direction_vector.x * self.z) * math.sin(teta)

        z = (point0.z * (direction_vector.x * direction_vector.x + direction_vector.y * direction_vector.y) -
             direction_vector.z * (point0.x * direction_vector.x + point0.y * direction_vector.y - direction_vector.x * self.x
                                   - direction_vector.y * self.y - direction_vector.z * self.z)) * (1 - math.cos(teta)) + self.z * math.cos(teta) + \
            (-point0.y * direction_vector.x + point0.x * direction_vector.y - direction_vector.y * self.x + direction_vector.x * self.y) * math.sin(teta)

        return Point(x, y, z)

    def getNormalPoint(self, direction_vector, point_vector, point):
        len_proection = Vector.scalar_product(direction_vector, point_vector)
        x = point.x + len_proection * direction_vector.x
        y = point.y + len_proection * direction_vector.y
        z = point.z + len_proection * direction_vector.z
        return Point(x, y, z)

class Vector(object):
    def __init__(self, point0, point1):
        import math
        self.x = point1.x - point0.x
        self.y = point1.y - point0.y
        self.z = point1.z - point0.z
        self.len = math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    @classmethod
    def scalar_product(cls, vector0, vector1):
        return vector0.x * vector1.x + vector0.y * vector0.z * vector1.z

    @classmethod
    def vector_product(cls, vector0, vector1):
        x = vector0.y * vector1.z - vector0.z * vector1.y
        y = vector0.z * vector1.x - vector0.x * vector1.z
        z = vector0.x * vector1.y - vector0.y * vector1.x
        point0 = Point(0, 0, 0)
        point1 = Point(x, y, z)
        return cls(point0, point1)

    def unit_vector(self):
        return Vector(Point(0, 0, 0), Point(self.x / self.len, self.y / self.len, self.z / self.len))

def getFiRiSi(Mi1, Oi, direction_vector, last_point):
    fi = Vector(last_point, Oi)
    ri = Vector(Oi, Mi1)
    si = Vector.vector_product(ri.unit_vector(), direction_vector)
    return fi, ri, si

def getNewConfiguration(coordinates, types, last_point):
    import math
    curr_coordinates = coordinates[:]
    for iAtom0 in range(0, len(coordinates) - 1, 3):
        iAtom1 = iAtom0 + 1
        iAtom2 = iAtom0 + 2

        if not ((types[iAtom0] == 'N' and types[iAtom1] == 'C' and types[iAtom2] == 'C') or
                    (types[iAtom0] == 'C' and types[iAtom1] == 'C' and types[iAtom2] == 'N')):
            continue

        M01 = Point(coordinates[iAtom0][0], coordinates[iAtom0][1], coordinates[iAtom0][2])
        M11 = Point(coordinates[iAtom1][0], coordinates[iAtom1][1], coordinates[iAtom1][2])
        direction_vector = Vector(M01, M11).unit_vector()

        for iAtom3 in range(iAtom2, len(coordinates), 1):
            M31 = Point(coordinates[iAtom3][0], coordinates[iAtom3][1], coordinates[iAtom3][2])
            O3 = M01.getNormalPoint(direction_vector, Vector(M01, M31), M01)
            f3, r3, s3 = getFiRiSi(M31, O3, direction_vector, last_point)
            teta = math.atan2(Vector.scalar_product(f3, s3) * r3.len, Vector.scalar_product(f3, r3) * r3.len)
            curr_point = M31.rotate(M31, direction_vector, teta)
            curr_coordinates[iAtom3] = (curr_point.x, curr_point.y, curr_point.z)

    return curr_coordinates

main()