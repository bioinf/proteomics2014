--------------------------
DESCRIPTION
--------------------------

Cyclic Coordinate Descent

Input:
PDB file with ATOM coordinates, coordinates of the target point (x, y, z)

Output:
PDF file with changed ATOM coordinates

--------------------------
USAGE
--------------------------

Example:  java -jar lab2.jar  1mbs.pdb  0.0 1.0 2.0

Usage: java -jar lab2.jar <in.pdb> x y z
    <in.pdb>               - input PDB file with ATOM coordinates
    x, y, z                - coordinates of the target point

Options:
    -h, --help, ?           Show this help
    -o <file>               Output file
