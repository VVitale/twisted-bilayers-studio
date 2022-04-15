# Twisted Bilayers Studio
This GitHub repository hosts the MATLAB codes developed at Imperial College for the calculation of the structural and electronic properties of twisted homo- and hetero- multilayer of different transition metal dichalcogenides.

## Authors
* Valerio Vitale (MATLAB codes for the generation of inital structures)
* Valerio Vitale and Kemal Atalar (MATLAB codes for tight-binding calculations of twisted bilayer TMDs)

## ASG - Atomic Structures Generator

Author
------

Valerio Vitale, Imperial College London, UK

* Email: vvitale@ic.ac.uk

Description
-----------

This module is a suite of MATLAB codes to generate the atomic structures of twisted homo- and hetero-bilayers of 2D
materials with a hexagonal unit cell, e.g. graphene, hBN, TMDs.
The input parameters can be specified in the `main.m` MATLAB file, both for homo and hetero bilayers. 
Several files can be generated:
1. `<rootname>.xyz` or `<rootname>.xsf` file for visualisation purposes
2. `positions.<rootname>.dat` file to be used as input (`positions.dat`) file for `TB_free.x`
3. `potential.<rootname>.dat` file to be used as input (`potential.dat`) file for `TB_free.x`
4. `lammps_positions.<rootname>.dat` file to be used as input file for `LAMMPS`

Usage
-----

You must have a MATLAB installed on your machine.

1. In MATLAB
* Change directory working directory to the one containing the `main.m` file relative to the system of interest (e.g. homobilayer)
* Simply type `main` in MATLAB Command Window
* Your output files can be found in your working directory

2. From terminal
* Open a new terminal window
* Change directory to the one containing the `main.m` relative to the system of interest (e.g. homobilayer)
* Type the following command:
	`$> matlab -nodisplay -nodesktop -nojvm -r "main ;exit"`
* Your output files can be found in your working directory

Notes
-----

For twisted homo-bilayers, e.g. twisted bilayer graphene (TBLG), the twist angle is specified by a pair of integers (n,m). 
The commensurate Moirè lattice is constructed from n,m and the lattice parameter of the single layer. 
The resulting Moirè cell is also hexagonal. A description of the input file is given in the header of the `main.m` file

For twisted hetero-bilayers, e.g. graphene on hBN, one has to specify a target angle and the ratio of the two lattice parameters. 
A description of the input parameters is given in the header of the `main.m` file



## TB - Tight-Binding

