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
NOTE: Tested on OpenFOAM 2.4.0, 3.0.0, 3.0.1, 4.0, 4.1 and 5.0 with gcc/Intel compilers.

Thanks to the contribution from @sanguinariojoe, the installation is also simplified:
```bash
of240 # change the version code according your OF installation
cd dugksFoam/src

# suggested by @sanguinariojoe
cp -r $WM_DIR/rules/$WM_ARCH$WM_COMPILER "$WM_DIR/rules/$WM_ARCH"MPI
export WM_COMPILER=MPI
sed -i "s/$WM_CXX/mpicxx/" $WM_DIR/rules/$WM_ARCH$WM_COMPILER/c++

./Allwmake
```

NOTE: the installation and usages in the `dugksFoam.pdf` in the `doc` directory is obsolete. The documentation will be rewritten.

## References
* [1] Z. Guo, K. Xu, R. Wang, Discrete unified gas kinetic scheme for all Knudsen number flows: low-speed isothermal case, [Phys. Rev. E, 88 (2013) 033305](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.88.033305).
* [2] Z. Guo, R. Wang, K. Xu, Discrete unified gas kinetic scheme for all Knudsen number flows. II. Thermal compressible case, [Phys. Rev. E, 91(2015) 033313.](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.033313)
* [3] L. Zhu, Z. Guo, K. Xu, Discrete unified gas kinetic scheme on unstructured meshes, [Comp. Fluids, 127(2016) 211-225](http://www.sciencedirect.com/science/article/pii/S0045793016000177)
* [4] L. Zhu, S. Chen, Z. Guo, dugksFoam: An open source OpenFOAM solver for the Boltzmann model equation, [Comp. Phys. Commun., 213(2017) 155-164](http://www.sciencedirect.com/science/article/pii/S0010465516303642)
* [5] L. Zhu, Z. Guo, Numerical study of nonequilibrium gas flow in a microchannel with a ratchet surface, [Phys. Rev. E. 95, (2017), 023113](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.95.023113)
