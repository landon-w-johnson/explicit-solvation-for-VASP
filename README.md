# Requirements

1. Python with NumPy





# How to Use This

1. Put your unsolvated `POSCAR` file in your working directory (whichever directory has `addExplicitSolvent.py`).

1. Put a `POSCAR` of your desired solvent molecule in your working directory with the name `SOLVENT`. Some samples are included in the `solvents` directory.

1. Run `python addExplicitSolvent.py`.

1. While `addExplicitSolvent.py` is running, you will be prompted to enter the density of your solvent in grams per mililiter. All other outputs are for debugging purposes.

1. After `addExplicitSolvent.py` is done running, you will have a new file called `POSCAR_WITH_SOLVENT` that contains your original molecule surrounded by randomly oriented solvent molecules at approximately your desired density. This file is ready for VASP.





# Known Bugs

1. This currently only works if the simulation cell in your `POSCAR` file is a rectangular prism.

1. This does not work with the Cartesian coordinate scheme yet.

1. It is not yet checking for collisions when placing solvent molecules, so VASP will very likely crash if you try to use a high density. You may want to visualize the output to make sure you don't have any atoms sitting on top of each other.