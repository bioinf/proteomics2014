2nd proteomics hometask
=======================
*problem3_solution.py* - this script solves loop closure problem with the help of CCD algorithm (problem 3).

Input:
- 1st parameter - pdb with input ATOM positions (only CA, C, N atoms are taken into consideration, everything else is ignored, different chains are ignored too)
- 2nd parameter - coordinates of last atom in chain (target positions of last atom in chain)
- 3rd parameter - output pdb filename. all data is copied from 1st pdb file, except atom coordinates from N, CA, C atoms - these gets modified (other atom coordinates are not modified if they present).

Output: pdb with new atom coordinates

Language and Notes
------------------

Script *problem3_solution.py* is written in Python 2.7, however, this is not actual hometask solution - I plan to solve problem 2 and to provide it as my 2nd proteomics hometask solution (although I'm lazy enough to write it...).

It runs once from N-terminal to C-terminal backbone atom. If target coordinates and current coordinates are still far from each other, calls script once again, untill epsilon=0.002 is reached or max_iterations=1000 is exceeded.

Currently problem 3 solution is fixed, and ready for checking.
