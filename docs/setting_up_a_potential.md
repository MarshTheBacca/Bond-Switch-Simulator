# Setting Up a Potential

There is already a template potential for you to use ![here](/testing/Si_potential.in)
Please note that as set in `Si_script.in`, the units used for the potentials are in electronic units.

## Bonds
`bond_style` is a ![LAMMPS command](https://docs.lammps.org/bond_style.html)
There are several styles to choose from, eg, ![harmonic](https://docs.lammps.org/bond_harmonic.html), ![morse](https://docs.lammps.org/bond_morse.html) etc.
It sets up the equation LAMMPS uses to simulate __bond__ interactions
`bond_coeff` is a ![LAMMPS command](https://docs.lammps.org/bond_coeff.html) that you use to give numeric values to the parameters defined in the bond_style. So for ![harmonic](https://docs.lammps.org/bond_harmonic.html), you need to specify _K_ and _r0_.
Please note __the first number for bond_coeff refers to the numerical ID of the bond type as defined in the data file__.
So, for example, if your data file has a `bond labels` section, you use the numerical ID, not the alpha-numeric label.


## Angles
`angle_style` is a ¬[LAMMPS command](https://docs.lammps.org/angle_style.html)
There are several styles to choose from, eg, ![harmonic](https://docs.lammps.org/angle_harmonic.html), ![cosine](https://docs.lammps.org/angle_cosine.html) etc.
It sets up the equation LAMMPS uses to simulate __angle__ interactions
`angle_coeff` is a ![LAMMPS command](https://docs.lammps.org/angle_coeff.html) that you use to give numeric values to the parameters defined in the bond_style. So for ![harmonic](https://docs.lammps.org/angle_harmonic.html), you need to specify _K_ and _theta0_.
Please note __the first number for angle_coeff refers to the numerical ID of the angle type as defined in the data file__.
So, for example, if your data file has a `angle labels` section, you use the numerical ID, not the alpha-numeric label.

## Advanced Potentials
If you wish to define your own custom potentials, please refer to the ![LAMMPS documentation](https://docs.lammps.org/)