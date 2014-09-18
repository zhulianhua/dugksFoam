/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fvDVM

Description

    Finite Volume Discrete Vecloty  Method. Solves the Boltzmann equation in 
    a disrec
    directions in a participating media, not including scatter.

    i.e. dictionary
    \verbatim
        fvDVMparas
        {
            xi_max        1000.0;          // Max discrete velocity
            xi_min       -1000.0;          // Min discrete velocity
            nDV               41;          // Number of discrete velocity, shoud be 4*Z + 1
        }

        gasProperties
        {
            R            80.0;             // Specific gas constant
            omega        0.7;              // VHS viscosity ~ Temperature index
            Tref         Tref [0 0 0 1 0 0 0] 275.0; // Reference temperature
            muRef        muRef [1 -1 -1 0 0 0 0] 1.0e-3;  // Reference viscosity
            Pr           0.75;             // Prantl number
        }
    \endverbatim

SourceFiles
    fvDVM.C

\*---------------------------------------------------------------------------*/

#ifndef fvDVM_H
#define fvDVM_H

#include "discreteVelocity.H"
#include "dimensionedScalar.H"
#include "calculatedMaxwellFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class discreteVelocity;
/*---------------------------------------------------------------------------*\
                           Class fvDVM Declaration
\*---------------------------------------------------------------------------*/

class fvDVM
:
    public IOdictionary
{
    // Private data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;

        //- Reference to the density field
        volScalarField& rhoVol_;

        //- Reference to the Macro velocity field
        volVectorField& Uvol_;

        //- Reference to the temperature field
        volScalarField& Tvol_;

        //- DVM model parameters dictionary
        const dictionary fvDVMparas_;

        //- Gas properties dictionary
        const dictionary gasProperties_;

        //- Number of discrete velocities in each direction
        const label nXiPerDim_;

        //- Total number of discrete velocities 
        label nXi_;

        label nXiX_;

        label nXiY_;

        label nXiZ_;

        //- Max discrete velocities
        const dimensionedScalar xiMax_;

        //- Min discrete velocities
        const dimensionedScalar xiMin_;

        //- grid space  in velocity space
        const dimensionedScalar dXi_;

        //- interg... value is  dXi^D, but dim is dXi^3
        const dimensionedScalar dXiCellSize_;

        //- use or not use macro flux
        const word macroFlux_;

        //- res
        const scalar res_;
        //- nStep for res
        const label checkSteps_;

        //- Gas: R
        const dimensionedScalar R_;

        //- Gas: omega
        const scalar omega_;

        //- Gas: Tref
        const dimensionedScalar Tref_;

        //- Gas: muRef
        const dimensionedScalar muRef_;

        //- Gas: Pr
        const scalar Pr_;

        //- Gas molecular innter degree of freedom
        const label KInner_;

        //- List of pointers to discrete velocities
        PtrList<discreteVelocity> DV_;

        //- surface density
        surfaceScalarField rhoSurf_;

        //- surface temperature
        surfaceScalarField Tsurf_;

        //- surface velocity
        surfaceVectorField Usurf_;

        //- surface heat flux vector
        surfaceVectorField qSurf_;

        //- surface stress tensor $\int \bm \xi_ \bm \xi_ d \bm \xi$
        surfaceTensorField stressSurf_;

        //- volume heat flux vector
        volVectorField qVol_;

        //- tau for cell
        volScalarField tauVol_;

        //- tau for face
        surfaceScalarField tauSurf_;


    // Private Member Functions
    //
        //- set discrete velocity
        void setDVgrid
        (
            scalarField& weights,
            scalarField& Xis,
            scalar xiMin,
            scalar xiMax,
            label nXi
        );

        //- Initialise
        //  Create the discrete velocity(DV) set and Initialize 
        //  each of the DV the with their discrete velocity value and weight
        void initialiseDV();

        //- set the outGoingByRho Field of calculatedMaxwellWall type rho Field
        //  do it only at the beginnig
        void setCalculatedMaxwellRhoBC();

        //- 1. Update h/g bar plus at cell center, per DV
        void updateGHbarPvol();

        //- 2.1 Update h/g bar at cell face (Interpolation), per DV
        void updateGHbarSurf();

        //- 2.2 Update the rho Boudndary Field, global DV
        //  Simply call the rhoVol_'s correctBoundaryConditions()
        void updateMaxwellWallRho(); 

        //- 2.3 Update the incoming maxwell wall DF, perDV
        void updateGHbarSurfMaxwellWallIn();

        //- 3. update Macros at cell faces, global DV
        void updateMacroSurf();

        //- 4. Update h/g at cell face h/g bar to h/g, per DV
        void updateGHsurf();

        //- 5. Update h/g tilde at cell center, per DV
        void updateGHtildeVol();

        //- 6. update Macros at cell centers, global DV
        void updateMacroVol();


        //- Disallow default bitwise copy construct
        fvDVM(const fvDVM&);

        //- Disallow default bitwise assignment
        void operator=(const fvDVM&);

    public:

        //- Runtime type information
        TypeName("fvDVM");


        // Constructors

        //- Construct from components
        fvDVM
        (
            volScalarField& rho, 
            volVectorField& U,
            volScalarField& T
        );

        //- Destructor
        virtual ~fvDVM();


        // Member functions

        // Edit

        //- Solve DVB equation(s), call the updateXxx functions.
        void evolution();

        // Access

        //- Co num, max and mean
        void getCoNum(scalar& maxCoNum, scalar& meanCoNum);

        //- discreteVelocity object for discreteVelocity I
        inline const discreteVelocity& DVi(const label I) const;

        //- discreteVelocity object for discreteVelocity xv, yv, zv
        inline const discreteVelocity& DVxyz
        (
            const label ix,
            const label iy,
            const label iz
        ) const;

        //- Density field at cell center
        inline volScalarField& rhoVol();

        //- Temperature field at cell center
        inline const volScalarField& Tvol() const;

        //- Macro velocity field at cell center
        inline const volVectorField& Uvol() const; 

        //- Heat flux field at cell center
        inline const volVectorField& qVol() const;

        //- tau at cell
        inline const volScalarField& tauVol() const;

        //- Density field at cell face
        inline const surfaceScalarField& rhoSurf() const;

        //- Temperature field at cell face
        inline const surfaceScalarField& Tsurf() const;

        //- Macro velocity field at cell face
        inline const surfaceVectorField& Usurf() const; 

        //- Heat flux field at cell face
        inline const surfaceVectorField& qSurf() const;
        //- stress tensor at cell face
        inline const surfaceTensorField& stressSurf() const;

        //- tau at cell face
        inline const surfaceScalarField& tauSurf() const;

        //- Total number of discrete velocities
        inline label nXi() const;

        //- Number of discrete velocities in each direction
        inline label nXiPerDim() const;

        inline dimensionedScalar xiMax() const;

        inline dimensionedScalar xiMin() const;

        inline dimensionedScalar R() const;

        inline scalar res() const;

        inline label checkSteps() const;

        inline scalar omega() const;

        inline dimensionedScalar Tref() const;

        inline dimensionedScalar muRef() const;

        inline scalar Pr() const;

        inline label KInner() const;

        //- aux. update tau from new T at vol/surf, call by updateMacroXxx
        template<template<class> class PatchType, class GeoMesh> 
        tmp<GeometricField<scalar, PatchType, GeoMesh> >
        updateTau
        (
            const GeometricField<scalar, PatchType, GeoMesh>& T, 
            const GeometricField<scalar, PatchType, GeoMesh>& rho
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvDVMI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
