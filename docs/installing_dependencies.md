# Installing Dependencies

## Packages

Install package dependencies with the following command:

```
sudo apt install build-essential clang-20 cmake ffmpeg git libfftw3-dev libjpeg-dev libomp-dev libopenmpi-dev libpng-dev libspdlog-dev ninja-build pkg-config
```

Link Clang 20 to Clang

```
ln -s /usr/bin/clang-20 /usr/bin/clang && ln -s /usr/bin/clang++-20 /usr/bin/clang++
```

## LAMMPS

CLone the [LAMMPS](https://www.lammps.org/) repository using `git clone https://github.com/lammps/lammps`

Enter the _lammps_ directory and execute the following commands:

```
mkdir build && cd build
cmake -C ../cmake/presets/basic.cmake -DBUILD_SHARED_LIBS=ON -DPKG_MOLECULE=yes -DPKG_EXTRA-MOLECULE=yes -DWITH_JPEG=yes -DWITH_PNG=yes -DWITH_FFMPEG=yes -DWITH_GZIP=yes  -DBUILD_MPI=yes -DBUILD_OMP=yes ../cmake
cmake --build . -j `nproc`
cmake --install .
```

By default, things will be installed to _~./local/lib_, and you need to add another line to _~/.bashrc_ `export CMAKE_PREFIX_PATH="$HOME/.local"`
This directory is searched by the CMake commands `find_library` and `find_path` when locating LAMMPS
