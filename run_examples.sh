#!/bin/bash

python3 main.py tests/li-ec/li-ec-01.inp -psf li-ec.psf -pdb li-ec.pdb
python3 main.py tests/fsi-tinker/fsi.inp -psf fsi-tinker.psf -pdb fsi-tinker.pdb -t
python3 main.py tests/fsi-tinker-drude/fsi.inp -psf fsi-npol.psf -pdb fsi-npol.pdb -t
python3 main.py tests/fsi-tinker-drude/fsi.inp -psf fsi-pol.psf -pdb fsi-pol.pdb -t -d
