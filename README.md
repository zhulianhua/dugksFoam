# dugksFoam
An OpenFOAM solver for Boltzmann model equation using discrete unified gas kinetic scheme [1,2,3].

## Key features
* Solving discrete velocity Boltzmann equation;
* Based on Shakhov-BGK collision model;
* Using discrete unified gas kinetic scheme (Asymptotic preserving property);
* A standard OpenFOAM solver;
* Arbitrary unstructured meshes;
* MPI parallel computing capability;
* 1D & 2D & 3D in a single solver;
* Post processing tools of OpenFOAM are ready to use;
* Various boundary condition types.

## Installation and documentation
```bash
of240 # change the version code according your OF installation
cd dugksFoam/src

# make own copy of wmake so we can override the 'CC' defination to mpicxx to support the velocity space decomposition MPI parallel computing ability.
# note: this workround may only works at OpenFOAM-2.4.0

cp `which wmake` .
sed -i '316s/makeType/makeType CC=mpicxx/' wmake
./Allwmake
```

## Test only on the following version of OpenFOAM  and platform
* OpenFOAM-2.4.0 with intel compilers
* OpenFOAM-2.4.0 with gnu compilers

See more details in `dugksFoam.pdf` in the `doc` directory, or download it [here](https://github.com/zhulianhua/dugksFoam/raw/master/doc/dugksFoam.pdf).

## References
* [1] Z.L. Guo, K. Xu, R.J. Wang, Discrete unified gas kinetic scheme for all Knudsen number flows: low-speed isothermal case, [Phys. Rev. E, 88 (2013) 033305](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.88.033305).
* [2] Z.L. Guo, R.J. Wang, K. Xu, Discrete unified gas kinetic scheme for all Knudsen number flows. II. Thermal compressible case, [Phys. Rev. E, 91(2015) 033313.](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.033313)
* [3] L.H. Zhu, Z.L. Guo, K. Xu, Discrete unified gas kinetic scheme on unstructured meshes, [Comp. Fluids, 127(2016) 211-225](http://www.sciencedirect.com/science/article/pii/S0045793016000177)
* [4] L.H. Zhu, S.Z. Chen, Z.L. Guo, dugksFoam: An open source OpenFOAM solver for the Boltzmann model equation, [Comp. Phys. Commun., 213(2017) 155-164](http://www.sciencedirect.com/science/article/pii/S0010465516303642)
