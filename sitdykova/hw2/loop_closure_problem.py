import sys
from Bio.PDB import *
import pdb_utils
from math import sqrt

eps = 0.001
max_iteration = 50

def dist(a, b):
    sum_of_squares = 0
    for i in range(0, 3):
        sum_of_squares += (a[i] - b[i]) * (a[i] - b[i])
    return sqrt(sum_of_squares)

def is_zero(x):
    return x < 0.0005

def get_rotation_axis(point1, point2, cur, target):
    v1 = cur - point1
    v2 = point2 - point1
    theta = Vector.normalized(v2)
    o = point1 + (v2.left_multiply((v1 * v2) / (v2 * v2)))
    ortho = cur - o
    ri = Vector.norm(ortho)
    if is_zero(ri):
        axis = (target - o) ** theta
        if is_zero(Vector.norm(axis)) or is_zero(Vector.norm((cur - o) ** theta)):
            return Vector(0, 0, 0)
        return Vector.normalized(axis)

    s = (Vector.normalized(ortho)) ** theta
    fi = target - o
    if not is_zero(abs(fi * s)):
        return theta
    return s


def rotate(theta, point1, c, cos_t, sin_t):
    v1 = c - point1
    o = point1 + (theta.left_multiply(v1 * theta))
    ortho = c - o
    ri = Vector.norm(ortho)
    if is_zero(ri):
        return c
    r = Vector.normalized(ortho)
    s = r ** theta
    res = o + s.left_multiply(sin_t * ri) + ortho.left_multiply(cos_t)
    return res

def get_S_coeffecients(theta, point1, cur, target):
    v1 = cur - point1
    o = point1 + (theta.left_multiply(v1 * theta))
    ri = cur - o
    fi = target - o
    si = Vector.normalized(ri) ** theta

    if is_zero(Vector.norm(ri)):
        return (0, 0, 0)

    a = (ri * ri) + (fi * fi)
    b = (ri * fi) * 2
    c = (si * fi) * 2 * sqrt(ri * ri)
    return Vector(a, b, c)

def CCD(atoms, target):
    iteration = 1
    distance = dist(atoms[-1], target)
    while iteration < max_iteration and distance > eps:
        print("Iteration " + str(iteration))
        for i in range(0, len(atoms) - 1):
            cur = Vector.copy(atoms[-1]);

            if dist(cur, target) < eps:
                break

            theta = get_rotation_axis(atoms[i], atoms[i + 1], cur, target)

            if is_zero(Vector.norm(theta)):
                continue

            s = get_S_coeffecients(theta, atoms[i], cur, target)

            if is_zero(s[1] ** 2 + s[2] ** 2):
                continue

            cos_t = s[1] / sqrt(s[1] ** 2 + s[2] ** 2)
            sin_t = s[2] / sqrt(s[1] ** 2 + s[2] ** 2)

            for j in range(i + 1, len(atoms)):
                atoms[j] = rotate(theta, atoms[i], atoms[j], cos_t, sin_t)

        iteration += 1
        distance = dist(atoms[-1], target)
        print("Distance to target: "+ str(distance))

if len(sys.argv) < 5:
    print("Usage: python loop_closure_problem.py <pdb file> x y z")
    exit()
else:
    pdb_file = sys.argv[1]
    target = Vector(sys.argv[2], sys.argv[3], sys.argv[4])
    print("Read pbd file...")
    structure, atoms = pdb_utils.get_vectors_from_pdb(pdb_file)
    print("CCD...")
    CCD(atoms, target)
    print("Write result to file out_" + pdb_file + "...")
    pdb_utils.write_pdb(structure, atoms, "out_" + pdb_file)
