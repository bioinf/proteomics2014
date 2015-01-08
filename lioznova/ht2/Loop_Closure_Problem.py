from Bio.PDB import *
import numpy
import sys
import os
import argparse
import timeit
from time import strftime

def get_rotation(prev_prev_coord, prev_coord, cur_coord, F):
	def get_vectors(axis, Mi_dash, F):
		ri = vector_to_axis(axis, Mi_dash)
		Oi = Mi_dash - ri
		if ri.norm() < 10**(-10): # Mi_dash lie on axis
			return
		ri_hat = ri.normalized()
		fi = F - Oi
		theta_hat = axis.normalized()
		si_hat = ri_hat ** theta_hat
		return (ri, fi, ri_hat, si_hat)

	def get_rotaxis(cos, sin, axis):
		vector = axis.copy()
		vector.normalize()
		(x, y, z) = vector.get_array()
		rotation_matrix = numpy.zeros((3, 3))
		rotation_matrix[0, 0] = (1 - cos) * x * x + cos
		rotation_matrix[0, 1] = (1 - cos) * x * y - sin * z
		rotation_matrix[0, 2] = (1 - cos) * x * z + sin * y
		rotation_matrix[1, 0] = (1 - cos) * x * y + sin * z
		rotation_matrix[1, 1] = (1 - cos) * y * y + cos
		rotation_matrix[1, 2] = (1 - cos) * y * z - sin * x
		rotation_matrix[2, 0] = (1 - cos) * x * z - sin * y
		rotation_matrix[2, 1] = (1 - cos) * y * z + sin * x
		rotation_matrix[2, 2] = (1 - cos) * z * z + cos
		return rotation_matrix

	prev_coord_v = prev_coord - prev_prev_coord
	cur_coord_v = cur_coord - prev_prev_coord
	F_v = F - prev_prev_coord
	axis = prev_coord_v
	Mi_dash = cur_coord_v
	vectors = get_vectors(axis, Mi_dash, F_v)
	flag = True
	if vectors:
		(ri, fi, ri_hat, si_hat) = vectors
	else:
		axis = Mi_dash ** F_v
		axis.normalize()
		(ri, fi, ri_hat, si_hat) = get_vectors(axis, Mi_dash, F_v)
		flag = False
	if fi * si_hat == 0: # planar case: fi * (ri ** theta) == 0
		axis = si_hat
		(ri, fi, ri_hat, si_hat) = get_vectors(axis, Mi_dash, F_v)
		flag = False
	b = 2 * ri.norm() * (fi * ri_hat)
	c = 2 * ri.norm() * (fi * si_hat)
	cos_alpha = b / numpy.sqrt(b * b + c * c)
	sin_alpha = c / numpy.sqrt(b * b + c * c)
	return (get_rotaxis(cos_alpha, sin_alpha, axis), flag)

def rotate(input_file_name, output_file_name, F, MAX_ITERATION_NUM, ACCURACY):
	structure = (PDBParser()).get_structure('DATA', input_file_name)
	atom_list = [a for a in structure.get_atoms()]

	first_point_index = None
	second_point_index = None
	start_rotation_index = None
	last_point_index = None
	for i in xrange(len(atom_list)):
		cur_atom = atom_list[i]
		if cur_atom.get_name() in ['CA', 'C', 'N']:
			if first_point_index == None:
				first_point_index = i
				continue
			if second_point_index == None:
				second_point_index = i
				start_rotation_index = i
				continue
			last_point_index = i

	iteration_num = 0
	while (iteration_num < MAX_ITERATION_NUM) and (((atom_list[last_point_index]).get_vector() - F).norm() > ACCURACY):
		prev_coord = atom_list[first_point_index].get_vector()
		prev_coord_index = first_point_index
		for i in xrange(second_point_index, last_point_index + 1):
			cur_atom = atom_list[i]
			cur_coord = cur_atom.get_vector()
			if cur_atom.get_name() in ['CA', 'C', 'N']:
				prev_prev_coord = prev_coord
				prev_prev_coord_index = prev_coord_index
				prev_coord = cur_coord
				prev_coord_index = i
				(rotation, flag) = get_rotation(prev_prev_coord, prev_coord, (atom_list[last_point_index]).get_vector(), F)
				if not flag:
					start_rotation_index = prev_prev_coord_index
				else:
					start_rotation_index = prev_coord_index
				translation = prev_prev_coord.get_array()
				for j in xrange(start_rotation_index + 1, len(atom_list)):
					v = atom_list[j].get_vector() - prev_prev_coord
					atom_list[j].set_coord(v.get_array())
					atom_list[j].transform(rotation, translation)
		iteration_num += 1

	io = PDBIO()
	io.set_structure(structure)
	io.save(output_file_name)

if __name__ == "__main__":
	start = timeit.timeit()
	if len(sys.argv) == 1:
		print "Usage:", sys.argv[0], "-i <name of pdb file> -o <output file name> -F <x y z of final point>"
		print "Example:", sys.argv[0], "-i in.pdb -o out.pdb -F 5.0 -10.0 0.0"
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='Cyclic Coordinate Descent for Loop Closure Problem.')
	parser.add_argument("-i", "--in_file", help="name of input pdb file", required=True)
	parser.add_argument("-o", "--out_file", help="name of output file", required=True)
	parser.add_argument('-F', '--coords', help="x y z", nargs='+', type=float, required=True)
	parser.add_argument('--accuracy', help = "accuracy of rotation", type=float, required=False, default=0.001)
	parser.add_argument('--max_iteration', help = "maximum number of iterations", type=int, required=False, default=100)

	args = parser.parse_args()
	input_file_name = args.in_file
	output_file_name = args.out_file
	if not os.path.isfile(input_file_name):
		print >> sys.stderr, "Not a file\t" + input_file
		exit(1)
	if not len(args.coords) == 3:
		print >> sys.stderr, "x y z expected, not\t" + args.coord
		exit(1)
	F = Vector(args.coords)
	MAX_ITERATION_NUM = args.max_iteration
	ACCURACY = args.accuracy
	rotate(input_file_name, output_file_name, F, MAX_ITERATION_NUM, ACCURACY)
	print "Ready! Time elapsed: ", timeit.timeit() - start, "sec"



