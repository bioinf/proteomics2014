import sys
from math import sqrt
#returns None if not atom, else returns tuple with coordinates
def parse_line(line):
    if len(line) < 54 or not line.startswith("ATOM  ") or line[12:16].strip() not in ["C", "CA", "N"]:
        return None
    return map(float, [line[30:38], line[38:46], line[46:54]])
#method reads data from file
def read_coordinates(pdb_file):
    pdb = open(pdb_file)
    result = filter(None, [ parse_line(line) for line in pdb])
    pdb.close()
    return result
def write_coordinates(pdb_input_file, pdb_output_file, modified_coordinates):
    pdb_output = open(pdb_output_file, 'w')
    with open(pdb_input_file) as pdb_input:
        for line in pdb_input:
            if parse_line(line) == None:
                pdb_output.write(line)
            else:
                current_coordinates = modified_coordinates.pop(0)
                pdb_output.write(
                    line[ : 30] +
                    "".join(map("{0:>8.3f}".format, current_coordinates)) +
                    line[54 : ]
                )
    pdb_output.close()
        #
#helper methods for vector manipulation
def dot(v1, v2): return sum(map(lambda(x, y): x * y, zip(v1, v2)))
def cross(v1, v2): return [v1[1] * v2[2] - v1[2] * v2[1], v1[2] * v2[0] - v1[0] * v2[2], v1[0] * v2[1] - v1[1] * v2[0]]
def vect(p1, p2): return map(lambda(x, y): x - y, zip(p2, p1))
def mult(v1, scalar): return map(lambda x: 1.0 * x * scalar, v1)
def plus(v1, v2): return map(lambda (x, y): x + y, zip(v1, v2))
#method rotates given point c relative to given axis (p1, p2) to given angle defined as (cost, sint)
def rotate(p1, p2, c, cost, sint):
    v1 = vect(p1, c)
    v2 = vect(p1, p2)
    theta = mult(v2, 1.0/sqrt(dot(v2, v2)))
    o = plus(p1, mult(v2, dot(v1, v2)/dot(v2, v2)))
    ortho = vect(o, c)
    ri = sqrt(dot(ortho, ortho))
    if ri == 0: return c
    r = mult(ortho, 1.0/ri)
    s = cross(r, theta)
    return plus(o, plus(mult(s, sint * ri), mult(ortho, cost)))
#method finds coefficients in minimized function for given target and current coordinates
def compute_ri_fi_si(p1, p2, current, target):
    v1 = vect(p1, current)
    v2 = vect(p1, p2)
    theta = mult(v2, 1.0/sqrt(dot(v2, v2)))
    o = plus(p1, mult(v2, dot(v1, v2)/dot(v2, v2)))
    ortho = vect(o, current)
    ri = sqrt(dot(ortho, ortho))
    r = (mult(ortho, 1.0/ri) if ri != 0 else ortho)
    s = cross(r, theta)
    return (ortho, vect(o, target), s)
#method finds coefficients in minimized function for given set of target_coordinates
def compute_s_coeffs(p1, p2, current_coordinates, target_coordinates):
    vectors = [
                compute_ri_fi_si(p1, p2, current_coordinates[i], target_coordinates[i])
                for i in range(len(current_coordinates))
              ]
    return reduce(
        lambda (a, b, c), (x, y, z): (a + x, b + y, c + z),
        [   (
                    dot(ri, ri) + dot(fi, fi),
                    dot(ri, fi) * 2,
                    dot(si, fi) * 2* sqrt(dot(ri, ri))
                ) for (ri, fi, si) in vectors
                if dot(ri, ri) != 0.0
            ], (0,0,0)
    )
def comp_all(all_coordinates, target_coordinates):
    for i in range(len(all_coordinates) - 2):
        num_of_coords = len(target_coordinates)
        current_coordinates = all_coordinates[- num_of_coords : ]
        (a, b, c) = compute_s_coeffs(all_coordinates[i], all_coordinates[i + 1], current_coordinates, target_coordinates)
        if b * b + c * c == 0:
            continue
        cost, sint = b/sqrt(b * b + c * c), c/sqrt(b * b + c * c)
        for j in range(i + 2, len(all_coordinates)):
            before = all_coordinates[j]
            all_coordinates[j] = rotate(all_coordinates[i], all_coordinates[i + 1], all_coordinates[j], cost, sint)
    return all_coordinates
if len(sys.argv) < 4:
    print "  call format: script.py <input_pdb> [[0.1,2.0,2.0],[2.3,4.0,5]] <output_pdb>\n    second parameter - array of coordinates of last n points (n>=1)"
    exit(1)
target_coordinates = eval(sys.argv[2])
result = comp_all(read_coordinates(sys.argv[1]), target_coordinates)
write_coordinates(sys.argv[1], sys.argv[3], result)
print("done")
