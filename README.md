libSCOF – Singularity-Constrained Octahedral Fields for Hexahedral Meshing
======

`libSCOF` is an implementation of [SCOF](https://www.meshing.rwth-aachen.de/publication/017/) \[[Liu et al. 2018](https://dl.acm.org/citation.cfm?id=3197517.3201344)\] distributed under GPLv3.

Note: this is a preliminary version without including the branching as in our paper. The complete version will be committed soon.

If you make use of `libSCOF` in your scientific work, please cite our paper. For your convenience,
you can use the following bibtex snippet:

    @article{Liu:2018:SOF:3197517.3201344,
    author = {Liu, Heng and Zhang, Paul and Chien, Edward and Solomon, Justin and Bommes, David},
    title = {Singularity-constrained Octahedral Fields for Hexahedral Meshing},
    journal = {ACM Trans. Graph.},
    issue_date = {August 2018},
    volume = {37},
    number = {4},
    month = jul,
    year = {2018},
    issn = {0730-0301},
    pages = {93:1--93:17},
    articleno = {93},
    numpages = {17},
    url = {http://doi.acm.org/10.1145/3197517.3201344},
    doi = {10.1145/3197517.3201344},
    acmid = {3201344},
    publisher = {ACM},
    address = {New York, NY, USA},
    keywords = {hexahedral meshing, integer-grid maps, octahedral fields, singularity graph},
   }

## What is SCOF?

SCOF is an algorithm for generating octahedral fields with the prescribed hexmeshable singularity graphs.

The input to the algorithm is a tetrahedral mesh together with singularity graph constraints in OpenVolumeMesh file format. For more details about the input file format, please refer to the Readme.txt in [InputFiles.] (https://www.meshing.rwth-aachen.de/publication/017/)

The output of the algorithm are a topological octahedral field encoded as matchings and field alignment as well as a smooth geometrical octahedral field in the form of rotation matrix.

## Compiling

libSCOF can be compiled independently resulting in a command line tool or compiled together with OpenFlipper. For more details about compiling, please see BUILDING.

## License

`libSCOF` is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version. See [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).


