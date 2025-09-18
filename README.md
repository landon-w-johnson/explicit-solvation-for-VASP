# Requirements

1. Python with NumPy





# How to Use This

1. Put your unsolvated `POSCAR` file in your working directory (whichever directory has `addExplicitSolvent.py`).

1. Put a `POSCAR` of your desired solvent molecule in your working directory with the name `SOLVENT` (watch out for Windows automatically adding the `.txt` extension). Some samples are included in the `solvents` directory.

1. Run `python addExplicitSolvent.py`.

1. While `addExplicitSolvent.py` is running, you will be prompted to enter the density of your solvent in grams per mililiter. All other outputs are for debugging purposes.

1. After `addExplicitSolvent.py` is done running, you will have a new file called `POSCAR_WITH_SOLVENT` that contains your original molecule surrounded by randomly oriented solvent molecules at approximately your desired density. This file is ready for VASP once you rename it to `POSCAR` in your VASP working directory.

1. Append the appropriate elemental species to your `POTCAR` based on the order that they appear in `POSCAR_WITH_SOLVENT`.





# Known Bugs

1. This will not work if your solvent molecule spans across the periodic boundaries in the file `SOLVENT`.

1. Your `POSCAR` and `SOLVENT` files cannot have any blank lines in the header (anything above the ion positions).

1. Your `POSCAR` file must have a blank line after the ion positions.





# Bug Reporting/Feature Implementation

1. Please contact me at landon.w.johnson@ndsu.edu if you notice any bugs or would like to have me implement certain features.

1. Let me know if there's a particular solvent molecule you want to use that isn't in the `solvents` subdirectory. I can optimize the structure of that molecule and add it to this repository.

1. This is currently set to prevent placing solvent atoms within 2 Angstroms of atoms within another molecule/material. This can be easily changed for user-defined cutoff distances.
