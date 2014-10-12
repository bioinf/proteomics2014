tree aligner hometask
=====================
progressive alignment tool.

Input: one or multiple .fasta files with amino acid sequences

Output: fasta file with gaps (now not actual fasta, string descriptions are temporarily lost)

Language, dependencies
----------------------

Project is written in Julia language (http://julialang.org/).

The main script depends on `ArgParse` package.
First, open interactive Julia console, by typing 'julia', then print command (it will install required package):
```
Pkg.add("AddParse")
```

Call example
------------
(relative to directory containing this README.md file)
```
julia src/main.jl -s data/BLOSUM62 -c "N" tests/sequences1.faa tests/res1.faa
```
There should be help message with arguments and their descriptions, `-s` key corresponds to scoring matrix file, `-c` - 1st letter of clustering type (`N`eighbour Joining, `U`PGMA, `WPGMA`).

First file should contain input protein strings, second file is there for writing down all results.

Current bottlenecks, future plans
---------------------------------

- score function with scorematrix seems to be slow

- I didn't tested in on big datasets

- I didn't follow all Performance Tips from Julia official docs.
