# psf-pdb-builder-python

Python script creating PSF and PDB files for running molecular dynamics simulations in NAMD, especially with usage of Drude particles.

# Requirements

Only standard library for Python 3.6 has been used, there are no further requirements.

# Running

```python3 main.py path_to_packmol_file [parameters]```

Available parameters:
* ```-psf file_name``` - produce PSF file with a given name
* ```-pdb file_name``` - produce PDB file with a given name
* ```-t``` or ```--tinker``` - for every residue .conn files are also provided
* ```-d``` or ```--drude``` - include declared Drude particles

Example running commands are presented in ```run_examples.sh```.

# Input files

To correctly run the program there are some input files needed:
* Packmol format input file which name is specified when running the program
* .xyz file with atoms position for the whole system
* .dat and .xyz files describing every residue of the system
* .conn files for every residue if ```-t``` option is set on

# Packmol input

Example:
```
seed 1030
tolerance 2.0
output li-ec-01.xyz
filetype xyz
structure li.xyz
  number 4
  inside cube -22.0 -22.0 -22.0 44.0
end structure
structure ec.xyz
  number 46
  inside cube -22.0 -22.0 -22.0 44.0
end structure
structure tfsi.xyz
  number 4
  inside cube -22.0 -22.0 -22.0 44.0
end structure
```

Only lines beginning with ```output```, ```structure``` and ```number``` are parsed.
* ```output``` - specifies name of the .xyz file containing atom positions in the whole system
* ```structure``` - for every residue of the system a section beginning with ```structure``` and ending with ```end``` should be specified.
In the line with ```structure``` keyword name of the .xyz file containing atom positions of given residue should be given.
Additionally ```number``` gives the number of repetitions of a given residue in the whole system

The order of ```structure``` declarations should be the same as they order in the .xyz file!

# .dat input file

For every residue declared as ```structure``` in the Packmol input file .dat file with the same name as specified .xyz file should be present. Example:
```
EC
O3 -0.53 15.9994
Ccc -0.05 12.0107
O3 -0.53 15.9994
C8 -0.05 12.0107
C8 -0.05 12.0107
O2 -0.53 15.9994
H1 0.13 1.00794
H1 0.13 1.00794
H1 0.13 1.00794
H1 0.13 1.00794
BONDS:
1 2
1 5
2 3
2 6
3 4
4 8
4 9
4 5
5 7
5 10
END
```

File structure description:
* first line - residue name (max 3 characters)
* for every atom one line should be specified in format ```atom_symbol_for_NAMD charge mass``` in the same order as they appear in the .xyz file
* at the end of the file ```BONDS``` section could be specified containing in separate lines indices of atoms forming bonds (the first atom has an index equal to 1), this section must be ended with ```END``` keyword
* if ```BONDS``` section is not specified the program will try to detect bonds automatically (if atoms are closer than 1.5 angstroms) or if program is launched with ```-t``` option bond information will be read from .conn files

To declare Drude atoms one needs to give ```d``` symbol and polarizability at the end of the line of the atom which should have the Drude atom, for example:
```
FSI
Sf 1.02 32.065 d 1.00
Ff -0.13 18.9984 d 1.00
Of -0.53 15.9994 d 1.00
Of -0.53 15.9994 d 1.00
Nf -0.66 14.0067
Sf 1.02 32.065 d 1.00
Ff -0.13 18.9984 d 1.00
Of -0.53 15.9994 d 1.00
Of -0.53 15.9994 d 1.00
```

In addition, when running program for given file but without specifying ```-d``` option, information about Drude particles will be skipped.

# .xyz file

For every residue declared in Packmol file a .xyz file containing list of atoms in it should be specified in standard .xyz format, example:
```
10

  O         -6.311309       -5.532314       -0.465320       
  C         -7.260238       -5.417400        0.432049       
  O         -7.345047       -4.170005        0.997326       
  C         -6.302795       -3.379001        0.350445       
  C         -5.446716       -4.423542       -0.336933       
  O         -7.870834       -6.396057        0.818066       
  H         -5.065063       -4.206253       -1.341663       
  H         -6.751815       -2.767502       -0.441364       
  H         -5.703414       -2.825454        1.136486       
  H         -4.593200       -4.729328        0.253094       
```

* first line - atoms number
* second line - for comments, skipped during parsing
* every next line - atom symbol and its Cartesian coordinates

# Tinker file (.conn)

If ```-t``` option was set on, for every residue there should exist a .conn file (with the same name as the .xyz file for this residue) in
the format of an output of Tinker analyze module, example is provided in [this file](tests/fsi-tinker/fsi.conn). From this file information
about bonds, angles and torsional angles is read.

# Examples description

Example input files are located in subdirectories in ```tests```, commands used to launch the program for them are specified in ```run_examples.sh``` file
* [li-ec](tests/li-ec) - basic example, without Drude particles and without Tinker output files with bonds specified in .dat files
* [fsi-tinker](tests/fsi-tinker) - example of system with Tinker output file
* [fsi-tinker-drude](tests/fsi-tinker-drude) - the same as the above, but with additionally declared Drude particles
