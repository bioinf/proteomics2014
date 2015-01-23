import sys
import os
import copy
import argparse

def read_fasta(fp):
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))

def get_seqs(input_file_name):
	text = {}
	with open(input_file_name) as fp:
		for name, seq in read_fasta(fp):
			text[name] = [seq]
	fp.close()
	return text

def read_matrix(fp):
	scoring_matrix = {}
	proteins = None
	with open(matrix_file_name) as fp:
		for line in fp:
			if (len(line) < 2) or (line[0] == "#"):
				continue
			if not proteins:
				proteins = (line.strip()).split()
				continue
			cur_line_char = (line.strip()).split()[0]
			for i in xrange(1, len((line.strip()).split())):
				w = int((line.strip()).split()[i])
				if not scoring_matrix.has_key(cur_line_char):
					scoring_matrix[cur_line_char] = {}
				scoring_matrix[cur_line_char][proteins[i-1]] = w
	return (scoring_matrix, proteins)

def score(arr1, arr2):
	result = 0
	for a1 in xrange(len(arr1)):
		for a2 in xrange(len(arr2)):
			result += scoring_matrix[proteins[a1]][proteins[a2]] * arr1[a1] * arr2[a2]
	return result

def freq_mat(s):
	fr = []
	for i in xrange(len(proteins)):
		fr_ = [0] * len(s[0])
		fr.append(fr_)
	for i in xrange(len(s)):
		for j in xrange(len(s[i])):
			fr[proteins.index(s[i][j])][j] += 1
	return fr


def align(p1_, p2_, distance_only):
	p1 = freq_mat(p1_)
	string_num1 = 0
	for i in xrange(len(proteins)):
		string_num1 += p1[i][0]
		p1[i].insert(0, 0)
	p1[len(proteins) - 1][0] = string_num1
	p2 = freq_mat(p2_)
	string_num2 = 0
	for i in xrange(len(proteins)):
		string_num2 += p2[i][0]
		p2[i].insert(0, 0)
	p2[len(proteins) - 1][0] = string_num2
	def space(num):
		arr = [0] * len(proteins)
		arr[-1] = num
		return arr
	def column(mat, pos):
		arr = [0] * len(proteins)
		for i in xrange(len(proteins)):
			arr[i] = mat[i][pos]
		return arr
	p1_len = len(p1[0])
	p2_len = len(p2[0])
	d = []
	for i in xrange(p1_len):
		d_ = [None] * p2_len
		d.append(d_)
	d[0][0] = 0
	for i in xrange(1, p1_len):
		d[i][0] = d[i-1][0] + score(column(p1, i), space(string_num2))
	for j in xrange (1, p2_len):
		d[0][j] = d[0][j-1] + score(space(string_num1), column(p2, j))
	for j in xrange(1, p2_len):
		for i in xrange(1, p1_len):
			d[i][j] = max(d[i-1][j] + score(column(p1, i), space(string_num2)), d[i][j-1] + score(space(string_num1), column(p2, j)), d[i-1][j-1] + score(column(p1, i), column(p2, j)))
	edit_dist = d[-1][-1]
	if distance_only:
		return edit_dist
	s_merged = [''] * (string_num1 + string_num2)
	(i, j) = (p1_len-1, p2_len-1)
	while (i > 0) and (j > 0):
		if (d[i][j] == d[i-1][j-1] + score(column(p1, i), column(p2, j))):
			for k in xrange(string_num1):
				s_merged[k] += (p1_[k][i-1])
			for k in xrange(string_num2):
				s_merged[string_num1 + k] += (p2_[k][j-1])
			(i, j) = (i-1, j-1)
		elif (d[i][j] == d[i-1][j] + score(column(p1, i), space(string_num2))):
			for k in xrange(string_num1):
				s_merged[k] += (p1_[k][i-1])
			for k in xrange(string_num2):
				s_merged[string_num1 + k] += ('*')
			(i, j) = (i-1, j)
		elif (d[i][j] == d[i][j-1] + score(space(string_num1), column(p2, j))):
			for k in xrange(string_num1):
				s_merged[k] += ('*')
			for k in xrange(string_num2):
				s_merged[string_num1 + k] += (p2_[k][j-1])
			(i, j) = (i, j-1)
	while (i > 0):
		for k in xrange(string_num1):
			s_merged[k] += (p1_[k][i-1])
		for k in xrange(string_num2):
			s_merged[string_num1 + k] += ('*')
		(i, j) = (i-1, j)
	while (j > 0):
		for k in xrange(string_num1):
			s_merged[k] += ('*')
		for k in xrange(string_num2):
			s_merged[string_num1 + k] += (p2_[k][j-1])
		(i, j) = (i, j-1)
	for i in xrange(len(s_merged)):
		s_merged[i] = s_merged[i][::-1]
	return s_merged

def initialize_dist(text):
	current_leaves = set(text.keys())
	profile_distances = {}
	for name in current_leaves:
		profile_distances[name] = {}
	for i in xrange(len(current_leaves)):
		for j in xrange(i+1, len(current_leaves)):
			p1_name = text.keys()[i]
			p2_name = text.keys()[j]
			p1_p2_dist = align(text[p1_name], text[p2_name], True)
			profile_distances[p1_name][p2_name] = p1_p2_dist
			profile_distances[p2_name][p1_name] = p1_p2_dist
	return (current_leaves, profile_distances)

def neighbor_joining(text_, profile_distances_, current_leaves_):
	text = copy.deepcopy(text_)
	profile_distances = copy.deepcopy(profile_distances_)
	current_leaves = copy.deepcopy(current_leaves_)
	while len(current_leaves) != 1:
		current_names = list(current_leaves)
		Q = {}
		for name in current_names:
			Q[name] = {}
		for i in xrange(len(current_names)):
			for j in xrange(i+1, len(current_names)):
				name1 = current_names[i]
				name2 = current_names[j]
				res = (len(current_names) - 2) * profile_distances[name1][name2]
				for k in xrange(len(current_names)):
					if k != i:
						res -= profile_distances[name1][current_names[k]]
					if k != j:
						res -= profile_distances[name2][current_names[k]]
				Q[name1][name2] = res
				Q[name2][name1] = res
		closest_dist = -100500
		closest_name1 = None
		closest_name2 = None
		for i in xrange(len(current_names)):
			for j in xrange(i+1, len(current_names)):
				if Q[current_names[i]][current_names[j]] > closest_dist:
					closest_dist = Q[current_names[i]][current_names[j]]
					closest_name1 = current_names[i]
					closest_name2 = current_names[j]
		new_profile_name = (closest_name1 + "\n" + closest_name2)
		new_profile = align(text[closest_name1], text[closest_name2], False)
		text[new_profile_name] = new_profile
		profile_distances[new_profile_name] = {}
		for cur_n in current_names:
			if (cur_n != closest_name1) and (cur_n != closest_name2):
				new_dist = 0.5 * (profile_distances[closest_name1][cur_n] + profile_distances[closest_name2][cur_n] - profile_distances[closest_name1][closest_name2])
				profile_distances[new_profile_name][cur_n] = new_dist
				profile_distances[cur_n][new_profile_name] = new_dist
		current_leaves.remove(closest_name1)
		current_leaves.remove(closest_name2)
		current_leaves.add(new_profile_name)
	final_name = current_leaves.pop()
	return (final_name, text[final_name])

def pgma(weighted, text_, profile_distances_, current_leaves_):
	text = copy.deepcopy(text_)
	profile_distances = copy.deepcopy(profile_distances_)
	current_leaves = copy.deepcopy(current_leaves_)
	while len(current_leaves) != 1:
		current_names = list(current_leaves)
		closest_dist = -100500
		closest_name1 = None
		closest_name2 = None
		for i in xrange(len(current_names)):
			for j in xrange(i+1, len(current_names)):
				if profile_distances[current_names[i]][current_names[j]] > closest_dist:
					closest_dist = profile_distances[current_names[i]][current_names[j]]
					closest_name1 = current_names[i]
					closest_name2 = current_names[j]
		new_profile_name = (closest_name1 + "\n" + closest_name2)
		new_profile = align(text[closest_name1], text[closest_name2], False)
		text[new_profile_name] = new_profile
		profile_distances[new_profile_name] = {}
		for cur_n in current_names:
			if (cur_n != closest_name1) and (cur_n != closest_name2):
				if weighted:
					new_dist = 0.5 * (profile_distances[closest_name1][cur_n] + profile_distances[closest_name2][cur_n])
				else:
					mod_name1 = len(text[closest_name1])
					mod_name2 = len(text[closest_name2])
					new_dist = float(mod_name1 * profile_distances[closest_name1][cur_n] + mod_name2 * profile_distances[closest_name2][cur_n]) / (mod_name1 + mod_name2)
				profile_distances[new_profile_name][cur_n] = new_dist
				profile_distances[cur_n][new_profile_name] = new_dist
		current_leaves.remove(closest_name1)
		current_leaves.remove(closest_name2)
		current_leaves.add(new_profile_name)
	final_name = current_leaves.pop()
	return (final_name, text[final_name])

def print_ans(names_, seqs, out_f):
	output_file = open(out_f, 'w')
	names = names_.split('\n')
	for i in xrange(len(seqs)):
		output_file.write(names[i] + "\n" + seqs[i] + "\n")
	output_file.close()

if __name__ == "__main__":
	if len(sys.argv) == 1:
		print "Usage:", sys.argv[0], "-o <output file> -i <input file> -s <scoring matrix> -m <NJ | UPGMA | WPGMA>"
		print "Please use the --help option to get more usage information."
		exit()

	parser = argparse.ArgumentParser(prog = sys.argv[0], description='msa')
	parser.add_argument("-i", "--in_file", help="name of input file", required=True)
	parser.add_argument("-o", "--out_file", help="name of output file", required=True)
	parser.add_argument("-s", "--score", help="scoring matrix file", required=True)
	parser.add_argument("-m", "--method", help="method to align", choices=["NJ", "UPGMA", "WPGMA"], required=True)

	args = parser.parse_args()
	input_file_name = args.in_file
	matrix_file_name = args.score
	output_file_name = args.out_file
	if not os.path.isfile(input_file_name):
		print >> sys.stderr, "Not a file\t" + input_file_name
		exit(1)
	if not os.path.isfile(matrix_file_name):
		print >> sys.stderr, "Not a file\t" + matrix_file_name
		exit(1)

	text = get_seqs(input_file_name)
	(scoring_matrix, proteins) = read_matrix(matrix_file_name)
	(current_leaves, profile_distances) = initialize_dist(text)

	if args.method == "NJ":
		(n, s) = neighbor_joining(text, profile_distances, current_leaves)
	elif args.method == "WPGMA":
		(n, s) = pgma(True, text, profile_distances, current_leaves) # wpgma
	elif args.method == "UPGMA":
		(n, s) = pgma(False, text, profile_distances, current_leaves) # upgma
	print_ans(n, s, output_file_name)

