# Installing Dependencies

## Open MPI
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

## LAMMPS

Download [LAMMPS](https://www.lammps.org/) from source, either as a tar.gz or from their [GitHub](https://github.com/lammps/lammps) using `git clone https://github.com/lammps/lammps`

Make a directory in the LAMMPS files called _build_

In the terminal, enter _build_ and execute the following commands:
```
cmake -C ../cmake/presets/basic.cmake -DBUILD_SHARED_LIBS=ON -DPKG_MOLECULE=yes -DPKG_EXTRA-MOLECULE=yes -DWITH_JPEG=yes -DWITH_PNG=yes -DWITHFFMPEG=yes -DWITH_GZIP=yes  -DBUILD_MPI=yes -DBUILD_OMP=yes ../cmake
cmake --build . -j X
cmake --install .
```
Where `X` is your number of cores (on Linux you can find this by executing `nproc`), this speeds up the build process.

By default, things will be installed to _~./local/lib_, and you need to add another line to _~/.bashrc_ `export PKG_CONFIG_PATH="$HOME/.local/lib/pkgconfig:$PKG_CONFIG_PATH"`
Now if you enter `pkg-config --list-all | grep liblammps` you should see liblammps is a recognised package.

## spdlog

Enter the following commands into the terminal:
```
git clone https://github.com/gabime/spdlog.git
cd spdlog && mkdir build && cd build
cmake .. && make -j
```
