#!/usr/bin/env python2

import sys
import os
from math import sqrt, pow
import traceback


class Atom():
    number = None
    x = None
    y = None
    z = None

    def __init__(self, splitted_line):
        
        self.number = int(splitted_line[1])
        self.x = float(splitted_line[6])
        self.y = float(splitted_line[7])
        self.z = float(splitted_line[8])
        self.coordinates = [self.x, self.y, self.z]

    def __iter__(self):
        for attr, value in self.__dict__.iteritems():
            yield attr, value

    def __repr__(self):
        return str(dict(self))

def cross_product(one, another):
        return ([one[1] * another[2] - one[2] * another[1],
                 one[2] * another[0] - one[0] * another[2],
                 one[0] * another[1] - one[1] * another[0]])

def plus(one, another):
        res = []
        for i in xrange(3):
            res.append(one[i] + another[i])
        return res

def get_vector(one, another):
        res = []
        for i in xrange(3):
            res.append(another[i] - one[i])
        return res

def length(one):
        res = 0.0
        for i in xrange(3):
            res += one[i] ** 2
        return sqrt(res)

def scalar_mult(one, scalar):
    res = []
    for i in range(3):
        res.append(one[i] * scalar)
    return res

def cdot(one, another):
        res = 0.0
        for i in xrange(3):
            res += (one[i] * another[i])
        return res

def normalize(one):
    return scalar_mult(one, 1.0 / length(one))

def parse_file(input_file):
    result_list_of_atoms = []
    with open(os.path.abspath(input_file)) as f:
        for line in f:
            if line.startswith("ATOM"):
                splitted_line = line.split()
                if splitted_line[2] in ["C", "CA", "N"]:
                    result_list_of_atoms.append(Atom(splitted_line).coordinates)
    return result_list_of_atoms

def calculate_o_point(point1, theta, vect1):
    return plus(point1, scalar_mult(theta, cdot(vect1, theta)))

def rotate_each_point(theta, point1, point2, sin, cos):
    vect1 = get_vector(point1, point2)
    o_point = calculate_o_point(point1, theta, vect1)
    orthogonal = get_vector(o_point, point2)
    r_i = length(orthogonal)
    if r_i < 10 ** (-10):
        return point2
    s = cross_product(normalize(orthogonal), theta)
    result = plus(o_point, plus(scalar_mult(s, sin * r_i), scalar_mult(orthogonal, cos)))
    return result

def get_s(theta, point1, current_coords, final_coords):
    vect1 = get_vector(point1, current_coords)
    o_point = calculate_o_point(point1, theta, vect1)
    r_i = get_vector(o_point, current_coords)
    f_i = get_vector(o_point, final_coords)
    s_i = cross_product(normalize(r_i) , theta)
    result = [0, 0, 0]
    if not cdot(r_i, r_i) == 0:
        result[0] += cdot(r_i, r_i) + cdot(f_i, f_i)
        result[1] += 2 * cdot(r_i, f_i)
        result[2] += 2 * sqrt(cdot(r_i, r_i)) * cdot(s_i, f_i)
    return result

def solve(input_file, output_file, final_coords_tuple, eps, iterations):
    coordinates = parse_file(input_file)
    iteration_counter = 0
    while True:
        iteration_counter += 1
        for i in range(len(coordinates) - 1):
            current_coord = coordinates[-1]
            if current_coord == list(final_coords_tuple):
                break
            theta = None
            vect1 = get_vector(coordinates[i], current_coord)
            vect2 = get_vector(coordinates[i], coordinates[i + 1])
            theta_new = normalize(vect2)
            o_point = plus(coordinates[i], scalar_mult(vect2, cdot(vect1, vect2) / cdot(vect2, vect2)))
            orthogonal = get_vector(o_point, current_coord)
            r_i = length(orthogonal)
            
            if r_i == 0:
                axis = cross_product(get_vector(o_point, final_coords_tuple), theta_new)
                if (length(axis) == 0 or length(cross_product(get_vector(o_point, current_coord), theta_new)) == 0):
                    theta = [0, 0, 0]
                else:
                    theta = normalize(axis)
            else:
                r = normalize(orthogonal)
                s = cross_product(r, theta_new)
                f_i = get_vector(o_point, final_coords_tuple)
                if not cdot(f_i, s) == 0:
                    theta = theta_new
                else:
                    theta = s
            
            if cdot(theta, theta) == 0:
                continue
            
            s = get_s(theta, coordinates[i], current_coord, final_coords_tuple)
            a = s[0]
            b = s[1]
            c = s[2]
            if b == 0 and c == 0:
                continue
            else:
                denominator = b * b + c * c
                cos = b / sqrt(denominator)
                sin = c / sqrt(denominator)
                for j in xrange(i + 1, len(coordinates)):
                    coordinates[j] = rotate_each_point(theta, coordinates[i], coordinates[j], sin, cos)
        dist_to_target = length(get_vector(coordinates[-1], final_coords_tuple))
        print "iteration: ", str(iteration_counter)
        print "distance: ", str(dist_to_target)
        if iteration_counter > iterations:
            print "Max number of iteration was reached"
            break
        if dist_to_target < eps:
            print "We are closer than epsilon to the solution"
            break
    i = 0
    with open(input_file, "r") as f1:
        with open(output_file, "w") as f2:
            for line in f1:
                if line.startswith("ATOM"):
                    splitted_line = line.split()
                    if splitted_line[2] in ["C", "CA", "N"]:
                        # can not invent new method of writing PDB because of shitty format specifications...
                        f2.write(line[:30] + '{0:>8.3f}'.format(coordinates[i][0]) + '{0:>8.3f}'.format(coordinates[i][1])
                                 + '{0:>8.3f}'.format(coordinates[i][2]) + line[54:])
                        i += 1
                    else:
                        f2.write(line.strip() + "\n")
                else:
                    f2.write(line.strip() + "\n")
                        



def main():
    help_string = """To execute this script, please, use our format of input command line.
\tpython2 loop_closure.py pdb.pdb out.pdb 0.0 1.0 2.0"""
    try:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        x_final = float(sys.argv[3])
        y_final = float(sys.argv[4])
        z_final = float(sys.argv[5])
        eps = 0.001
        iterations = 100
        solve(input_file, output_file, (x_final, y_final, z_final), eps, iterations)
    except :
        print traceback.format_exc()
        print help_string

if __name__ == "__main__":
    main()
