#!/usr/bin/python

import sys
import logging
from math import sqrt

from Bio.PDB import *

threshold = 0.001


def read_input(file):
    parser = PDBParser()
    structure, result = parser.get_structure('input', file), []
    for atom in structure.get_atoms():
        if atom.get_id() in ["C", "N", "CA"]:
            result.append(atom.get_vector())
    return structure, result


def get_rotation(point1, point2, cur, target):
    v1, v2 = cur - point1, point2 - point1
    theta = Vector.normalized(v2)
    o = point1 + (v2.left_multiply((v1 * v2) / (v2 * v2)))
    ortho = cur - o
    if Vector.norm(ortho) < threshold:
        temp = (target - o) ** theta
        return Vector(0, 0, 0) if (Vector.norm(temp) < threshold) or \
                                  ((Vector.norm((cur - o) ** theta)) < threshold) else Vector.normalized(temp)
    else:
        temp = (Vector.normalized(ortho)) ** theta
        return theta if not (abs((target - o) * temp)) < threshold else temp


def rotate(theta, point1, c, cos, sin):
    v1 = c - point1
    o = point1 + (theta.left_multiply(v1 * theta))
    ortho = c - o
    if Vector.norm(ortho) < threshold:
        return c
    else:
        s = Vector.normalized(ortho) ** theta
        res = s.left_multiply(sin * Vector.norm(ortho)) + ortho.left_multiply(cos) + o
        return res


def do_job(target, atoms):
    eps = 0.001
    max_iter = 50

    def dist(a, b):
        res = 0
        for i in range(3):
            res += (a[i] - b[i]) * (a[i] - b[i])
        return sqrt(res)

    def get_coeff(theta, point, cur, target):
        temp = point + (theta.left_multiply((cur - point) * theta))
        ri, fi = cur - temp, target - temp
        si = Vector.normalized(ri) ** theta
        return (0, 0, 0) if Vector.norm(ri) < threshold else Vector((ri * ri) + (fi * fi), (ri * fi) * 2, (si * fi) * 2 * sqrt(ri * ri))

    iter = 1
    distance = dist(atoms[-1], target)

    while iter < max_iter and distance > eps:
        logging.info("Iteration " + str(iter))

        for i in range(len(atoms) - 1):
            cur = Vector.copy(atoms[-1])

            if dist(cur, target) < eps:
                break

            theta = get_rotation(atoms[i], atoms[i + 1], cur, target)

            if Vector.norm(theta) < threshold:
                continue

            coef = get_coeff(theta, atoms[i], cur, target)
            if (coef[1] ** 2 + coef[2] ** 2) < threshold:
                continue

            for j in range(i + 1, len(atoms)):
                atoms[j] = rotate(theta, atoms[i], atoms[j], (coef[1] / sqrt(coef[1] ** 2 + coef[2] ** 2)),
                                  coef[2] / sqrt(coef[1] ** 2 + coef[2] ** 2))

        logging.info("Finish iteration")
        iter += 1
        distance = dist(atoms[-1], target)
        logging.info("Current distance " + str(distance))


def write_answer(file, structure, coords):
    i = 0

    for atom in structure.get_atoms():
        if atom.get_id() in ["C", "N", "CA"]:
            atom.set_coord(coords[i])
            i += 1

    io = PDBIO()
    io.set_structure(structure)
    io.save(file)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    if len(sys.argv) < 6:
        sys.stderr.write("USAGE: main.py <input_file> <x> <y> <z> <output_file>\n")
        sys.exit(1)

    in_file = sys.argv[1]
    target = Vector(sys.argv[2], sys.argv[3], sys.argv[4])
    out_file = sys.argv[5]

    logging.info("Read problem")
    structure, atoms = read_input(in_file)

    logging.info("Do job")
    do_job(target, atoms)

    logging.info("Write results")
    write_answer(out_file, structure, atoms)

    logging.info("Finish")
