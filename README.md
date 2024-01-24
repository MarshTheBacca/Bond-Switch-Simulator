# LAMMPS Network Monte Carlo

This is a fork of the Wilson Group's [2D Network Monte Carlo](https://github.com/WilsonGroupOxford/Network-Monte-Carlo) repository that integrates [LAMMPS](https://github.com/lammps/lammps) to allow for faster energy evaluation and other 2D structures:

* Silica Bilayers (Si2O3)
* Triangle Rafts (SiO2)
* Tersoff Graphene (C)
* Boron Nitride (BN)

## Installing Dependencies

### Open MPI
Download [Open MPI](https://www.open-mpi.org/)

Extract to a folder, make a new directory called _build_

In the terminal, enter _build_ and execute the following commands:
```
./configure 2>&1 | tee config.out
make -j 4 all 2>&1 | tee make.out
make install 2>&1 | tee install.out
```

You may have to install python dependencies if you encounter errors, like I had to install _sphinx_rtd_theme_ using `pip3 install sphinx_rtd_theme` (your command may be just `pip`, or if you're in an anaconda environment, `conda`)

If you still encounter an error with _sphinx_rtd_theme_, then enter "pip3 show sphinx_rtd_theme" to find where it is installed, then edit _~/.bashrc_ and add the line `export PYTHON_PATH="<path>:$PYTHON_PATH"`

Then execute `source ~/.bashrc` to update the environment variable (and it will automatically execute everytime you log on)

### LAMMPS

Download [LAMMPS](https://www.lammps.org/) from source, either as a tar.gz or from their [GitHub](https://github.com/lammps/lammps) using `git clone https://github.com/lammps/lammps`

Make a directory in the LAMMPS files called _build_

In the terminal, enter _build_ and execute the following commands:
```
cmake ../cmake -DBUILD_SHARED=ON
cmake --build . -j X --target install
```
Where `X` is your number of cores (on Linux you can find this by executing `nproc`), this speeds up the build process.

By default, things will be installed to _~./local/lib_, and you need to add another line to _~/.bashrc_ `export PKG_CONFIG_PATH="$HOME/.local/lib/pkgconfig:$PKG_CONFIG_PATH"`
Now if you enter `pkg-config --list-all | grep liblammps` you should see liblammps is a recognised package.

## Credit

Credit must be given to [Oliver Whitaker](https://github.com/oliwhitg) for implementing LAMMPS and [David Ormrod Morley](https://github.com/dormrod) for the origial NetMC code.
