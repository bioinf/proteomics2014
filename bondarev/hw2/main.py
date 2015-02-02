#!/usr/bin/env python
from __future__ import print_function, division

import sys
import logging
from itertools import islice
from math import sqrt

import numpy
import Bio.PDB as pdb

range = xrange


class LoopStructurePredictor:
    MAX_ITERATIONS_COUNT = 100

    def __init__(self, initial_structure):
        self.structure = initial_structure
        self.atoms = list(self.structure.get_atoms())
        self.backbone_indexes = [i for i, atom in enumerate(self.atoms) if atom.get_name() in ('CA', 'C', 'N')]

    def update_structure(self, target_coords):
        for iteration_index in range(self.MAX_ITERATIONS_COUNT):
            logging.info('Start iteration {}. Distance: {:.4e}'.format(iteration_index, self.distance(target_coords)))
            if self._is_close_enough(target_coords):
                break

            for prev_index, curr_index in self._paired_backbone_indices():
                rotation_matrix, start_rotation_index = self._get_rotation(prev_index, curr_index, target_coords)
                self._rotate_arm_end(self._coords(prev_index), rotation_matrix, start_rotation_index)
        else:
            logging.warning('Maximum number of iterations reached. '
                            'Result distance {:.4e}'.format(self.distance(target_coords)))
        return self.structure

    def distance(self, target_coords):
        return (self._loop_tail_coords() - target_coords).norm()

    def _is_close_enough(self, target_coords):
        return self.distance(target_coords) < 1e-6

    def __getitem__(self, atom_index):
        return self.atoms[atom_index]

    def _coords(self, atom_index):
        return self[atom_index].get_vector()

    def _loop_tail_coords(self):
        return self._coords(self.backbone_indexes[-1])

    def _paired_backbone_indices(self):
        # pair all atom indexes except last pair
        return zip(self.backbone_indexes, islice(self.backbone_indexes, 1, len(self.backbone_indexes) - 1))

    def _rotate_arm_end(self, prev_coord, rotation_matrix, start_rotation_index):
        translation = prev_coord.get_array()
        for pos in range(start_rotation_index + 1, len(self.atoms)):
            v = self._coords(pos) - prev_coord
            self[pos].set_coord(v.get_array())
            self[pos].transform(rotation_matrix, translation)

    def _get_rotation(self, prev_index, curr_index, target_coords):
        prev_coords = self._coords(prev_index)
        curr_coords = self._coords(curr_index)
        start_rotation_index = curr_index

        target_vector = target_coords - prev_coords
        rotation_axis = curr_coords - prev_coords
        m_point = self._loop_tail_coords() - prev_coords

        vectors = self._get_system_vectors(rotation_axis, m_point, target_vector)
        if vectors is not None:
            r, f, r_normd, s_normd = vectors
        else:
            rotation_axis = (m_point ** target_vector).normalized()
            r, f, r_normd, s_normd = self._get_system_vectors(rotation_axis, m_point, target_vector)
            start_rotation_index = prev_index

        if f * s_normd == 0:
            rotation_axis = s_normd
            r, f, r_normd, s_normd = self._get_system_vectors(rotation_axis, m_point, target_vector)
            start_rotation_index = prev_index

        b = 2 * r.norm() * (f * r_normd)
        c = 2 * r.norm() * (f * s_normd)
        cos_alpha = b / sqrt(b * b + c * c)
        sin_alpha = c / sqrt(b * b + c * c)
        return self._get_rotation_matrix(cos_alpha, sin_alpha, rotation_axis), start_rotation_index

    @staticmethod
    def _get_system_vectors(rotation_axis, m_point, target_vector):
        r = pdb.vector_to_axis(rotation_axis, m_point)  # the perpendicular projection m_point to rotation_axis
        o = m_point - r                                 # corresponded rotation axis vector

        if r.norm() < 1e-9:  # m_point on rotation axis
            return
        r_normd = r.normalized()
        f = target_vector - o
        theta_norm = rotation_axis.normalized()
        s_normd = r_normd ** theta_norm
        return r, f, r_normd, s_normd

    @staticmethod
    def _get_rotation_matrix(cos_a, sin_a, axis):
        x, y, z = axis.normalized().get_array()
        return numpy.array([
            [(1 - cos_a) * x * x + cos_a,     (1 - cos_a) * y * x - sin_a * z, (1 - cos_a) * z * x + sin_a * y],
            [(1 - cos_a) * x * y + sin_a * z, (1 - cos_a) * y * y + cos_a,     (1 - cos_a) * z * y - sin_a * x],
            [(1 - cos_a) * x * z - sin_a * y, (1 - cos_a) * y * z + sin_a * x, (1 - cos_a) * z * z + cos_a]
        ])


def usage():
    return ('Usage:\n'
            '  main.py <pdb> <fixed target coords: x y z> <output_filename>')


def main():
    logging.basicConfig(level=logging.INFO)
    if len(sys.argv) != 6:
        logging.error('Wrong count of arguments')
        print(usage())
        sys.exit()

    _, pdb_filename, x, y, z, out_filename = sys.argv

    fixed_target = pdb.Vector(x, y, z)
    logging.info('Read pdb file')
    structure = pdb.PDBParser().get_structure('DATA', pdb_filename)
    logging.info('Finish read pdb file')
    predictor = LoopStructurePredictor(structure)
    logging.info('Start calc result structure')
    predictor.update_structure(fixed_target)
    logging.info('Save result structure')
    io = pdb.PDBIO()
    io.set_structure(predictor.structure)
    io.save(out_filename)
    logging.info('Done.')


if __name__ == "__main__":
    main()