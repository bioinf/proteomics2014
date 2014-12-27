import sys
from math import sqrt
#returns None if not atom, else returns tuple with coordinates
def parse_line(line):
    if len(line) < 54 or not line.startswith("ATOM  ") or line[12:16].strip() not in ["C", "CA", "N"]:
        return None
    return map(float, [line[30:38], line[38:46], line[46:54]])

# method reads data from file
def read_coordinates(pdb_file):
    pdb = open(pdb_file, 'r')
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
def length(v): return sqrt(dot(v, v))
def unit(v):  return mult(v, 1.0/length(v))

#method rotates given point c relative to given axis (p1, p1+theta) to given angle defined as (cost, sint)
def rotate(theta, p1, c, cost, sint):
    v1 = vect(p1, c)
    o = plus(p1, mult(theta, dot(v1, theta)))
    ortho = vect(o, c)
    ri = length(ortho)
    if ri == 0: return c
    r = unit(ortho)
    s = cross(r, theta)
    return plus(o, plus(mult(s, sint * ri), mult(ortho, cost)))

def get_rotation_axis(p1, p2, current, target):
    v1 = vect(p1, current)
    v2 = vect(p1, p2)
    theta = unit(v2)
    o = plus(p1, mult(v2, dot(v1, v2)/dot(v2, v2)))
    ortho = vect(o, current)
    ri = length(ortho)
    if ri == 0:
        axis = cross(vect(o, target), theta)
        if (length(axis) == 0 or length(cross(vect(o, current), theta)) == 0):
            return [0, 0, 0]
        return unit(axis)
    r = unit(ortho)
    s = cross(r, theta)
    fi = vect(o, target)
    if dot(fi, s) != 0:
        return theta
    return s
    #there is another case, when fi and ri are in the same plane with theta.
    #at this point we take s as rotation axis


# method finds coefficients in minimized function for given target and current coordinates
# returns triplet of (theta, ri, fi, si) for given point p1
def compute_ri_fi_si(theta, p1, current, target):
    v1 = vect(p1, current)
    o = plus(p1, mult(theta, dot(v1, theta)))
    ri = vect(o, current)
    fi = vect(o, target)
    si = cross(unit(ri), theta)
    return (ri, fi, si)


# method finds coefficients in minimized function for given set of target_coordinates
def compute_s_coeffs(theta, p1, current_coordinate, target_coordinate):
    vectors = [compute_ri_fi_si(theta, p1, current_coordinate, target_coordinate)]

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
def comp_all(all_coordinates, target_coordinate, epsilon = 0.002, max_iterations = 1000):
    current_iteration = 0
    while (1):
        current_iteration += 1
        for i in range(len(all_coordinates) - 1):
            current_coordinate = all_coordinates[-1]
            if current_coordinate == target_coordinate:
                break
            theta = get_rotation_axis(all_coordinates[i], all_coordinates[i + 1], current_coordinate, target_coordinate)
            if dot(theta, theta) == 0:
                continue
            (a, b, c) = compute_s_coeffs(theta, all_coordinates[i], current_coordinate, target_coordinate)
            if b * b + c * c == 0:
                continue
            cost, sint = b/sqrt(b * b + c * c), c/sqrt(b * b + c * c)
            for j in range(i + 1, len(all_coordinates)):
                before = all_coordinates[j]
                all_coordinates[j] = rotate(theta, all_coordinates[i], all_coordinates[j], cost, sint)
        print(
            "N-terminal to C-terminal: iteration {0} passed,\n \
            difference in target coordinates and current coordinates for last atom is {1}\
            ".format(
                current_iteration,
                length(vect(all_coordinates[-1], target_coordinate))
            )
        )
        if (current_iteration >= max_iterations):
            print("Number of iterations exceeds max_iterations constant")
            break
        if (length(vect(all_coordinates[-1], target_coordinate)) <= epsilon):
            print("Number of iterations is less than epsilon={0}".format(epsilon))
            break
    return all_coordinates

def show_help():
    print("\'call me maybe' format: script.py <input_pdb> [0.1,2.0,2.0] <output_pdb>\n\
            \tsecond parameter - target coordinates of last point\n")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        show_help()
        exit(1)
    try:
        target_coordinate = eval(sys.argv[2])
        result = comp_all(read_coordinates(sys.argv[1]), target_coordinate)
    except SyntaxError:
        print("dude, check your syntax - all your rockets crushed and the sky are fallen")
        show_help()
        exit(1)
    except IOError:
        print("problem reading input file, check if exists")
        show_help()
        exit(1)
    except ValueError:
        print("got some problems with input pdb file")
        show_help()
        exit(1)
    else:
        write_coordinates(sys.argv[1], sys.argv[3], result)
        print("Good news, everyone! We've got a very special delivery today. \nIt's a brand new PDB file named " + sys.argv[3])
