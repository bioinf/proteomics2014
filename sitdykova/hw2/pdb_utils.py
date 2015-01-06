from Bio.PDB import *

def get_vectors_from_pdb(pdbfile):
    parser = PDBParser()
    header = parser.get_header
    structure = parser.get_structure('input', pdbfile)
    res = []
    for atom in structure.get_atoms():
        if atom.get_id() in ["C", "N", "CA"]:
            res.append(atom.get_vector())
    return structure, res

def write_pdb(structure, coords, new_file):
    i = 0
    for atom in structure.get_atoms():
        if atom.get_id() in ["C", "N", "CA"]:
            atom.set_coord(coords[i])
            i += 1
    io = PDBIO()
    io.set_structure(structure)
    io.save(new_file)