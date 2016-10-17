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
# change from g++/icpc to mpicxx to support velocity space decomposition MPI parallel computing ability.
cp $WM_DIR/rules/$WM_ARCH$WM_COMPILER/c++ $WM_DIR/rules/$WM_ARCH$WM_COMPILER/c++.bak # make a backup.
sed -i "s/$WM_CXX/mpicxx/"  $WM_DIR/rules/$WM_ARCH$WM_COMPILER/c++ 
./Allwmake
```

See more details in `dugksFoam.pdf` in the `doc` directory, or download it [here](https://github.com/zhulianhua/dugksFoam/raw/master/doc/dugksFoam.pdf).

## References
* [1] Z.L. Guo, K. Xu, R.J. Wang, Discrete unified gas kinetic scheme for all Knudsen number flows: low-speed isothermal case, [Phys. Rev. E, 88 (2013) 033305](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.88.033305).
* [2] Z.L. Guo, R.J. Wang, K. Xu, Discrete unified gas kinetic scheme for all Knudsen number flows. II. Thermal compressible case, [Phys. Rev. E, 91(2015) 033313.](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.033313)
* [3] L.H. Zhu, Z.L. Guo, K. Xu, Discrete unified gas kinetic scheme on unstructured meshes, arXiv preprint arXiv:1503.07374, (2015).
