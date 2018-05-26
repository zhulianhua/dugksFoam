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

\*---------------------------------------------------------------------------*/

#include "foam_defs.h"
#include "zeroBoundaryVectors.H"
#include "volFields.H"

#if FOAM_MAJOR <= 3
    #define BOUNDARY_FIELD_REF boundaryField()
#else
    #define BOUNDARY_FIELD_REF boundaryFieldRef()
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(zeroBoundaryVectors, 0);
}


// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

Foam::zeroBoundaryVectors::zeroBoundaryVectors(const fvMesh& mesh)
:
    MeshObject<fvMesh, Foam::MoveableMeshObject, zeroBoundaryVectors>(mesh),
    pVectors_
    (
        IOobject
        (
            "LeastSquaresP",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    ),
    nVectors_
    (
        IOobject
        (
            "LeastSquaresN",
            mesh_.pointsInstance(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh_,
        dimensionedVector("zero", dimless/dimLength, vector::zero)
    )
{
    calcLeastSquaresVectors();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::zeroBoundaryVectors::~zeroBoundaryVectors()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zeroBoundaryVectors::calcLeastSquaresVectors()
{
    if (debug)
    {
        Info<< "zeroBoundaryVectors::calcLeastSquaresVectors() :"
            << "Calculating least square gradient vectors"
            << endl;
    }

    const fvMesh& mesh = mesh_;

    // Set local references to mesh data
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceScalarField& w = mesh.weights();
    const surfaceScalarField& magSf = mesh.magSf();


    // Set up temporary storage for the dd tensor (before inversion)
    symmTensorField dd(mesh_.nCells(), symmTensor::zero);

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector d = C[nei] - C[own];
        symmTensor wdd = (magSf[facei]/magSqr(d))*sqr(d);

        dd[own] += (1 - w[facei])*wdd;
        dd[nei] += w[facei]*wdd;
    }

#if FOAM_MAJOR <= 3
    surfaceVectorField::GeometricBoundaryField& blsP =
        pVectors_.boundaryField();
#else
    surfaceVectorField::Boundary& blsP =
        pVectors_.BOUNDARY_FIELD_REF;
#endif

    forAll(blsP, patchi)
    {
        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.patch().faceCells();

        // Build the d-vectors
        vectorField pd(p.delta());

        if (pw.coupled())
        {
            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                dd[faceCells[patchFacei]] +=
                    ((1 - pw[patchFacei])*pMagSf[patchFacei]/magSqr(d))*sqr(d);
            }
        }
        else
        {
            // We do nothing here, since we ignore the boundary faces when 
            // calc the grad
            // This is how does this BC differs the least Square
            //forAll(pd, patchFacei)
            //{
                //const vector& d = pd[patchFacei];

                //dd[faceCells[patchFacei]] +=
                    //(pMagSf[patchFacei]/magSqr(d))*sqr(d);
            //}
        }
    }


    // Invert the dd tensor
    const symmTensorField invDd(inv(dd));


    // Revisit all faces and calculate the pVectors_ and nVectors_ vectors
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        vector d = C[nei] - C[own];
        scalar magSfByMagSqrd = magSf[facei]/magSqr(d);

        pVectors_[facei] = (1 - w[facei])*magSfByMagSqrd*(invDd[own] & d);
        nVectors_[facei] = -w[facei]*magSfByMagSqrd*(invDd[nei] & d);
    }

    forAll(blsP, patchi)
    {
        fvsPatchVectorField& patchLsP = blsP[patchi];

        const fvsPatchScalarField& pw = w.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        const fvPatch& p = pw.patch();
        const labelUList& faceCells = p.faceCells();

        // Build the d-vectors
        vectorField pd(p.delta());

        if (pw.coupled())
        {
            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                patchLsP[patchFacei] =
                    ((1 - pw[patchFacei])*pMagSf[patchFacei]/magSqr(d))
                   *(invDd[faceCells[patchFacei]] & d);
            }
        }
        else
        {
            forAll(pd, patchFacei)
            {
                const vector& d = pd[patchFacei];

                patchLsP[patchFacei] =
                    pMagSf[patchFacei]*(1.0/magSqr(d))
                   *(invDd[faceCells[patchFacei]] & d);
            }
        }
    }

    if (debug)
    {
        Info<< "zeroBoundaryVectors::calcLeastSquaresVectors() :"
            << "Finished calculating least square gradient vectors"
            << endl;
    }
}


bool Foam::zeroBoundaryVectors::movePoints()
{
    calcLeastSquaresVectors();
    return true;
}


// ************************************************************************* //
