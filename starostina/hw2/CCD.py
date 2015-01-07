import sys
from math import sqrt 

class Vector:
    def __init__(self, x = 0, y = 0, z = 0, p1 = None, p2 = None):
        self.x = x
        self.y = y
        self.z = z
        if p1 and p2:
            self.x = p2.x - p1.x
            self.y = p2.y - p1.y
            self.z = p2.z - p1.z

    def dot_product(self, v):
        return self.x * v.x + self.y * v.y + self.z * v.z

    def cross_product(self, v):
        return Vector(self.y * v.z - self.z * v.y, self.z * v.x - self.x * v.z, self.x * v.y - self.y * v.x)

    def num_product(self, num):
        return Vector(self.x * num, self.y * num, self.z * num)
    
    def __add__(self, v):
        return Vector(self.x + v.x, self.y + v.y, self.z + v.z)

    def length(self):
        return sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def unit(self):
        length = self.length()
        return Vector(self.x / length, self.y / length, self.z / length)

def getO(theta, p, c):
    v = Vector(p1 = p, p2 = c)
    return theta.num_product(v.dot_product(theta)) + p

class Atom:
    def __init__(self, x, y, z, ind = None, type = None):
        self.index = ind
        self.type = type
        self.x = x
        self.y = y
        self.z = z

    def rotate(self, theta, p, cos, sin):
        o = getO(theta, p, self)
        r = Vector(p1 = o, p2 = self)
        ortho_len = r.length()
        if ortho_len < 0.01:
            return
        r = r.unit()
        s = r.cross_product(theta)
        res = o + s.num_product(sin * ortho_len) + r.num_product(cos * ortho_len)
        self.x = res.x
        self.y = res.y
        self.z = res.z

def read_pdb(infile):
    atoms = []
    for line in open(infile):
        if not line.startswith('ATOM'):
            continue
        fields = line.split()
        atom = Atom(float(fields[6]), float(fields[7]), float(fields[8]), int(fields[1]), fields[2])
        atoms.append(atom)
    return atoms

def write_res_pdb(infile, outfile, atoms):
    fin = open(infile)
    fout = open(outfile, 'w')
    i = 0
    for line in fin:
        if not line.startswith('ATOM'):
            fout.write(line.strip() + '\n')
        else:
            fout.write(line[:30] + '{0:>8.3f}'.format(atoms[i].x) + '{0:>8.3f}'.format(atoms[i].y) + '{0:>8.3f}'.format(atoms[i].z) + line[54:])
            i += 1


def get_NC_indices(atoms):
    indices = []
    for i in range(len(atoms)):
        if atoms[i].type in ['CA', 'C', 'N']:
            indices.append(i)
    return indices

def get_rotation_axis(p1, p2):
    axis = Vector(p1 = p1, p2 = p2)
    return axis.unit()

def calc_ri_fi_si(theta, p, curr_end_pos, target_end_pos):
    o = getO(theta, p, curr_end_pos)
    ri = Vector(p1 = o, p2 = curr_end_pos)
    fi = Vector(p1 = o, p2 = target_end_pos)
    si = ri.unit().cross_product(theta)
    return (ri, fi, si)

def calc_S_coeffs(theta, p, curr_end_pos, target_end_pos):
    ri, fi, si = calc_ri_fi_si(theta, p, curr_end_pos, target_end_pos)
    a = ri.length() + fi.length()
    b = ri.dot_product(fi)
    c = si.dot_product(fi) * 2 * sqrt(ri.length())
    return (a, b, c)

def CCD(atoms, target_end_pos, accuracy, max_iteration):
    NC_indices = get_NC_indices(atoms)
    iteration = 0
    curr_dist = 100500
    while iteration <= max_iteration and curr_dist > accuracy:
        for i in range(len(NC_indices) - 1):
            curr_end_pos = atoms[NC_indices[-1]]
            p1 = atoms[NC_indices[i]]
            p2 = atoms[NC_indices[i + 1]]
            theta = get_rotation_axis(p1, p2)
            a, b, c = calc_S_coeffs(theta, p1, curr_end_pos, target_end_pos)
            if b == 0 and c == 0:
                continue
            cos = b / sqrt(b * b + c * c)
            sin = c / sqrt(b * b + c * c)

            for j in range(NC_indices[i] + 1, len(atoms)):
                atoms[j].rotate(theta, p1, cos, sin)

        iteration += 1
        curr_dist = Vector(p1 = atoms[NC_indices[-1]], p2 = target_end_pos).length()
        print 'Distance to target position after ' + str(iteration) + ' iteration is ' + str(curr_dist)
    return atoms

if len(sys.argv) != 6:
    print "Usage: ./CCD.py in.pdb out.pdb target_x target_y target_z"
    sys.exit(-1)

input_pdb = sys.argv[1]
output_pdb = sys.argv[2]
target_atom = Atom(float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
atoms = CCD(read_pdb(input_pdb), target_atom, 0.05, 100)
write_res_pdb(input_pdb, output_pdb, atoms)
