/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::discreteVelocity

Description
    Discrete velocity object

SourceFiles
    discreteVelocity.C

\*---------------------------------------------------------------------------*/

#ifndef discreteVelocity_H
#define discreteVelocity_H

#include "fvCFD.H"
#include "calculatedMaxwellFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvDVM;

/*---------------------------------------------------------------------------*\
                    Class discreteVelocity Declaration
\*---------------------------------------------------------------------------*/

class discreteVelocity
{
public:
        static const word discreteVelocityPrefix;

private:

    // Private data

        //- Refence to the owner fvDVM object
        fvDVM& dvm_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the time
        const Time& time_;

        //- weight
        const scalar weight_;

        //- discrete velocity 
        const dimensionedVector xi_;

        //- My DV Id
        const label myDVid_;

        //- 
        const label symXtargetDVid_;
        const label symYtargetDVid_;
        const label symZtargetDVid_;

        //- volume data
        volScalarField gTildeVol_;

        volScalarField hTildeVol_;

        volScalarField gBarPvol_;

        volScalarField hBarPvol_;


        //- surface data
        surfaceScalarField gSurf_;

        surfaceScalarField hSurf_;


        //- gradient of g/h bar B
        volVectorField gBarPgrad_;
        volVectorField hBarPgrad_;


        //- Global DV id - incremented in constructor
        static label DVid_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        discreteVelocity(const discreteVelocity&);

        //- Disallow default bitwise assignment
        void operator=(const discreteVelocity&);


public:

    // Constructors

        //- Construct form components
        discreteVelocity
        (
            fvDVM& dvm,
            const fvMesh& mesh,
            const Time& time,
            const scalar weight,
            const dimensionedVector xi,
            const label DVid,
            const label symXtargetDVid,
            const label symYtargetDVid,
            const label symZtargetDVid
        );

    //- Destructor
    ~discreteVelocity();


    // Member functions

        // Edit


            //- Initialise DF to equilibrium
            void initDFtoEq();

            //- set BC type for DF
            void setBCtype();

            void initBoundaryField();

            //- Update h/g bar plus at cell center
            void updateGHbarPvol();

            //- Update h/g bar at cell face (Interpolation)
            void updateGHbarSurf();

            //- Update the incoming maxwell type DF
            //  called after the execution of the updateMaxwellWallRho()
            void updateGHbarSurfMaxwellWallIn();

            //- Update the incomming type DF at symmetry plane
            void updateGHbarSurfSymmetryIn();

            //- Update h/g at cell face h/g bar to h/g
            void updateGHsurf();

            //- Update h/g tilde at cell center
            void updateGHtildeVol();


            //- get equilibrium using Maxwell equilibrium
            //  Used by initialization, once
            template <template<class> class PatchType, class GeoMesh> 
            Foam::tmp
            <
                Foam::GeometricField<Foam::scalar, PatchType, GeoMesh> 
            >
            equilibriumMaxwell
            (
                const GeometricField<scalar, PatchType, GeoMesh>& rho,
                const GeometricField<vector, PatchType, GeoMesh>& U,
                const GeometricField<scalar, PatchType, GeoMesh>& T
            );

            //- get gEq and hEq using Maxwell equilibrium
            //- Used every evolvetion step
            template <template<class> class PatchType, class GeoMesh> 
            void equilibriumShakhov
            (
                GeometricField<scalar, PatchType, GeoMesh>& gEq,
                GeometricField<scalar, PatchType, GeoMesh>& hEq,
                const GeometricField<scalar, PatchType, GeoMesh>& rho,
                const GeometricField<vector, PatchType, GeoMesh>& U,
                const GeometricField<scalar, PatchType, GeoMesh>& T,
                const GeometricField<vector, PatchType, GeoMesh>& q
            );

            //- get the equilibrium at boundary patches
            void equilibriumMaxwell
            (
                fvsPatchScalarField& geq,
                fvsPatchScalarField& heq,
                const fvPatchScalarField&  rho,
                const fvPatchVectorField&  U,
                const fvPatchScalarField&  T
            );

            //- get the equilibrium at each boundary face
            scalar equilibriumMaxwellByRho(const vector U, const scalar T);


        // Access

            //- Return xi
            inline const dimensionedVector xi() const;

            //- Return weight
            inline scalar weight() const;

            //- Return \tilde g 
            inline const volScalarField& gTildeVol() const;

            //- Return \tilde h
            inline const volScalarField& hTildeVol() const;

            //- Return \tilde g 
            inline const volScalarField& gBarPvol() const;

            //- Return \tilde h
            inline const volScalarField& hBarPvol() const;

            //- Return  g at faces
            inline const surfaceScalarField& gSurf() const;

            //- Return  h at faces
            inline const surfaceScalarField& hSurf() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "discreteVelocityI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
