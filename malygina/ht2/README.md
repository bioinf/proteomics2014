2nd proteomics hometask
=======================
*problem3_solution.py* - this script solves loop closure problem with the help of CCD algorithm (problem 3).

Input:
- 1st parameter - pdb with input ATOM positions (only CA, C, N atoms are taken into consideration, everything else is ignored, different chains are ignored too)
- 2nd parameter - array of coordinates of last atoms in chain (target positions of last n atoms, n>=1)
- 3rd parameter - output pdb filename. all data is copied from 1st pdb file, except atom coordinates from N, CA, C atoms - these gets modified.

Output: pdb with new atom coordinates

Language and Notes
------------------

Script *problem3_solution.py* is written in Python 2.7, however, this is not actual hometask solution - I plan to solve problem 2 and to provide it as my 2nd proteomics hometask solution.

Currently problem 3 solution has some unsolved issues, in particular, *situation when all atoms lay in the same plane isn't processed properly*.
