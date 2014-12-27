Падает на загрузке матрицы. Беру BLOSUM62.txt:

```
ERROR: BoundsError()
 in getindex at string.jl:625
 in readMatrix at /Users/pavel/DEV/biocad/other/proteomics2014/malygina/ht1/src/DataReader.jl:42
 in main at /Users/pavel/DEV/biocad/other/proteomics2014/malygina/ht1/src/main.jl:57
 in include at boot.jl:245
 in include_from_node1 at loading.jl:128
 in process_options at client.jl:285
 in _start at client.jl:354
while loading /Users/pavel/DEV/biocad/other/proteomics2014/malygina/ht1/src/main.jl, in expression starting on line 65
```
