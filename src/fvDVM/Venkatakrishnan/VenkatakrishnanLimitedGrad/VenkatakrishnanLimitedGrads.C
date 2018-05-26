/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
#include "VenkatakrishnanLimitedGrad.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "volFields.H"
#include "fixedValueFvPatchFields.H"

#if FOAM_MAJOR <= 3
    #define BOUNDARY_FIELD_REF boundaryField()
#else
    #define BOUNDARY_FIELD_REF boundaryFieldRef()
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if FOAM_MAJOR <= 2
namespace Foam
{
namespace fv
{
#endif

    makeFvGradScheme(VenkatakrishnanLimitedGrad)

#if FOAM_MAJOR <= 2
}
}
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
Foam::tmp<Foam::volVectorField>
Foam::fv::VenkatakrishnanLimitedGrad<Foam::scalar>::calcGrad
(
    const volScalarField& vsf,
    const word& name
) const
{
    const fvMesh& mesh = vsf.mesh();

    tmp<volVectorField> tGrad = basicGradScheme_().calcGrad(vsf, name);

    if (VenSlope_.k() < SMALL)
    {
        return tGrad;
    }

    volVectorField g = tGrad();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const volVectorField& C = mesh.C();
    const surfaceVectorField& Cf = mesh.Cf();
    const scalarField& V = mesh.V();

    scalarField maxVsf(vsf.internalField());
    scalarField minVsf(vsf.internalField());

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        scalar vsfOwn = vsf[own];
        scalar vsfNei = vsf[nei];

        maxVsf[own] = max(maxVsf[own], vsfNei);
        minVsf[own] = min(minVsf[own], vsfNei);

        maxVsf[nei] = max(maxVsf[nei], vsfOwn);
        minVsf[nei] = min(minVsf[nei], vsfOwn);
    }


#if FOAM_MAJOR <= 3
    const volScalarField::GeometricBoundaryField& bsf = vsf.boundaryField();
#else
    const volScalarField::Boundary bsf = vsf.boundaryField();
#endif

    forAll(bsf, patchi)
    {
        const fvPatchScalarField& psf = bsf[patchi];

        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();

        if (psf.coupled())
        {
            const scalarField psfNei(psf.patchNeighbourField());

            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                scalar vsfNei = psfNei[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
        else
        {
            forAll(pOwner, pFacei)
            {
                label own = pOwner[pFacei];
                scalar vsfNei = psf[pFacei];

                maxVsf[own] = max(maxVsf[own], vsfNei);
                minVsf[own] = min(minVsf[own], vsfNei);
            }
        }
    }

    maxVsf -= vsf;
    minVsf -= vsf;

    // create limiter
    volScalarField limiter(
        IOobject
        (
            "limiter",
            mesh.time().timeName(),
            mesh, 
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zeroes", dimensionSet(0,0,0,0,0,0,0), 0)
    );
    // scalarField limiter(vsf.internalField().size(), 1.0);

    scalar limiterTemp = 1.0;

    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // owner side
        limiterTemp = 
            VenSlope_.limitFace
            (
                V[own],
                maxVsf[own],
                minVsf[own],
                (Cf[facei] - C[own])&g[own]
            );

        limiter[own] = (limiter[own] < limiterTemp)?limiter[own]:limiterTemp;

        // neighbour side
        limiterTemp = 
            VenSlope_.limitFace
            (
                V[nei],
                maxVsf[nei],
                minVsf[nei],
                (Cf[facei] - C[nei])&g[nei]
            );
        limiter[nei] = limiter[nei] < limiterTemp?limiter[nei]:limiterTemp;
    }

    forAll(bsf, patchi)
    {
        const labelUList& pOwner = mesh.boundary()[patchi].faceCells();
        const vectorField& pCf = Cf.boundaryField()[patchi];

        forAll(pOwner, pFacei)
        {
            label own = pOwner[pFacei];

            limiterTemp = 
                VenSlope_.limitFace
                (
                    V[own],
                    maxVsf[own],
                    minVsf[own],
                    (pCf[pFacei] - C[own])&g[own]
                );
            limiter[own] = (limiter[own] < limiterTemp)?limiter[own]:limiterTemp;
        }
    }

    if (fv::debug)
    {
        Info<< "gradient limiter for: " << vsf.name()
            << " max = " << gMax(limiter)
            << " min = " << gMin(limiter)
            << " average: " << gAverage(limiter) << endl;
    }

    // g.internalField() *= limiter;
    g *= limiter;
    g.correctBoundaryConditions();
    gaussGrad<scalar>::correctBoundaryConditions(vsf, g);

    return tGrad;
}


// ************************************************************************* //

template<>
Foam::tmp<Foam::volTensorField>
Foam::fv::VenkatakrishnanLimitedGrad<Foam::vector>::calcGrad
(
    const volVectorField& vsf,
    const word& name
) const
{
    //NOTE by lhzhu, nothing done here for vector field 
    //const fvMesh& mesh = vsf.mesh();
    tmp<volTensorField> tGrad = basicGradScheme_().calcGrad(vsf, name);
    return tGrad;
}
