## Input file

All the parameters are read by ignoring all characters after the first white space. The lines the parameters appear on determine what internal variable they are assigned to, so you cannot change the format or order of the file, but you can change titles and dashed lines. The program has validation on the input data and will tell you if you've done something wrong.

__Fields marked with a '*' are not yet implemented__

| Name | Description | Validation |
| :--- | :--- | :--- |
| Output Folder | The directory to store output files relative to the directory containing netmc.x | String |
| Output File Prefix | The prefix each file will include, eg 'test' will produce 'test_A_crds.dat' | String |
| Input Folder | The directory to read input files from relative to the directory containing netmc.x | String |
| Input File Prefix | The prefix each input file includes | String |
| __Create Network From Scratch?*__ | Determins if the program reads an existing network or produces a new crystalline network | String 'true' or 'false' |
| __Number of Rings*__ | The number of rings to produce if creating a network from scratch | Integer > 1 |
| Min Ring Size | The maximum allowed ring size during the simulation excluding any 'fixed rings' | Integer > 2 |
| Max Ring Size | The minimum allowed ring size during the simulation excluding any 'fixed rings' | Integer >= Min Ring Size |
| Enable Fixed Rings? | Switches on the 'fixed rings' functionality of the program | String 'true' or 'false' |
| __Enable OpenMPI?*__ | Switches on the use of multi-thread processing to speed up simulation | String 'true' or 'false' |
| Structure Type | Determins what potential and LAMMPS data file is read, eg, 'Silicene' will read from si_potential.in and si.data | String 'Graphene' 'Silicene' __'TriangleRaft*'__ __'Bilayer*'__ __'BoronNitride*'__ |
| Random Seed | The seed used to generate random numbers in the program | Integer >= 0 |
| Selection Process | The method used to determine which bonds will be switched | String 'Random' 'Weighted' |
| Weighted Decay |  The exponential decay factor used if using a 'Weighted' selection process | Float |
| Thermalisation Temperature (10^x) | The temperature at which the system is thermalised at | Float |
| Annealing Start Temperature (10^x) | The temperature at which the annealing process starts at | Float |
| Annealing End Temperature (10^x) | The temperature at which the annealing process ends at | Float |
| Thermalisation Steps | The number of steps for the thermalisation process | Integer >= 0 |
| Annealing Steps | The number of steps for the annealing process | Integer >= 0 |
| Maximum Bond Length | The maxmimum bond length allowed for nodes involved in a switch move | Float > 0 |
| Maximum Bond Angle | The maximum bond angle allowed for nodes involved in a switch move | 0 < Float < 360 |
| Analysis Write Interval | How often analysis data will be written to all_stats.csv, if 0, no analysis written | Integer >= 0 |
| Write a Movie File? | Determines if movie.mpg is written, a video of the simulation happening (takes ~15x longer) | String 'true' or 'false' AND total steps <= 2000 |
