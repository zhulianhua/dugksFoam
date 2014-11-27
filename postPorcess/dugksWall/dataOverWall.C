/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Application
    dugksFoam

Description
    Discrete Unified Gas Kinetic Scheme(Zhaoli Guo, Kun Xu) solver.
    Lianhua Zhu

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvm.H"
#include "constants.H"
using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readTimeControls.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    // Read Fields
    surfaceVectorField qSurf
    (
        IOobject
        (
            "qSurf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    surfaceTensorField stressSurf
    (
        IOobject
        (
            "stressSurf",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Get wall pathcid
    label wallPatchID = mesh.boundaryMesh().findPatchID("cylinderWall");
    vectorField Cf = mesh.Cf().boundaryField()[wallPatchID];
    vectorField Sf = mesh.Sf().boundaryField()[wallPatchID];
    vectorField q = qSurf.boundaryField()[wallPatchID];
    tensorField stress = stressSurf.boundaryField()[wallPatchID];
    vectorField F(stress&Sf/mag(Sf));
    scalar ad = 0;
    forAll(Cf,facei)
    {
        scalar x = Cf[facei].component(0);
        scalar y = Cf[facei].component(1);
        scalar theta = 0;
        if(y>0)
        {

            if(x<0) ad=0;
            else ad = 180;
            theta = ad - 180.0/M_PI*Foam::atan(y/x);
            Info<<Cf[facei].component(0) << "\t" << Cf[facei].component(1);
            Info<<"\t"<< theta;
            Info<<"\t"<<(q[facei]&Sf[facei])/mag(Sf[facei]);
            Info<<"\t"<<(F[facei]&Sf[facei])/mag(Sf[facei]);
            Info<<"\t"<<magSqr(F[facei]- (F[facei]&Sf[facei])*Sf[facei]/magSqr(Sf[facei])) << endl;
        }
    }
    return 0;
}

// ************************************************************************* //
