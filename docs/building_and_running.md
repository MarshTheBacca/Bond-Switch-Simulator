# Building and Running

## Building

The program can be built using __cmake__ like so:

```
mkdir build; cd build
cmake ../ 
cmake --build .
```

You can use `-DEXEC_OUTPUT_PATH=<path>` in the second command to specify where you would like _bond_switch_simulator.exe_ to be located. You may want this to be in the _run_ folder (set path to _../run/_).

## Running

The program can be run using

```
./bond_switch_simulator.exe [-d]
```

The `-d` option is there to enable debugging messages. I wouldn't use this for lengthy simulations due to the ASCII pictograms logged for every single move in the simulation, which could potentially take up a large amount of storage.
