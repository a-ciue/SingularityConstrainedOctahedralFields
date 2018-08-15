libSCOF – Singularity-Constrained Octahedral Fields for Hexahedral Meshing
======

`libSCOF` is an implementation of [SCOF](https://www.meshing.rwth-aachen.de/publication/017/) \[[Liu et al. 2018] distributed under GPLv3.

Note: this is a preliminary version without including the branching as in our paper. The complete version will be committed soon.

If you make use of `libSCOF` in your scientific work, please cite our paper. 

## What is SCOF?

SCOF is an algorithm for generating octahedral fields with the prescribed hexmeshable singularity graphs.

The input to the algorithm is a tetrahedral mesh together with singularity graph constraints in OpenVolumeMesh file format. For more details about the input file format, please refer to the Readme.txt in [InputFiles.] () Additional constraints that are necessary to uniquely specify the field topology for models of higher genus, independent sub-singularity graphs or cavities, are not stored in the files yet.

The output of the algorithm are a topological octahedral field encoded as matchings and field alignment as well as a smooth geometrical octahedral field in the form of rotation matrix.

## Compiling

libSCOF can be compiled independently resulting in a command line tool or compiled together with OpenFlipper. For more details about compiling, please see BUILDING.

## License

`libSCOF` is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version. See [http://www.gnu.org/licenses/](http://www.gnu.org/licenses/).


