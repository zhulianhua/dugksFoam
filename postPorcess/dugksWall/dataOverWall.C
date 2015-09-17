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
#include <iomanip>
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
    scalar pm = 1;
    scalar theat = 0;
    scalar r = 0.01;

    Info << "x\ty\tq\tp\tsigma\n";
    std::cout << std::setprecision(8);
    std::cout << std::setw(12);
    std::cout << std::fixed;
    std::cout << std::scientific;
    forAll(Cf,facei)
    {
        scalar x = Cf[facei].component(0);
        scalar y = Cf[facei].component(1);
        scalar theta = 0;
        if(y>0)
        {
            if(y>2*r && x < 0) { ad =   0; pm = -1;}
            if(y>2*r && x > 0) { ad = 180; pm = -1;}
            if(y<2*r && x > 0) { ad = 180; pm = -1;}
            if(y<2*r && x < 0) { ad = 360; pm = -1;}

            //if(x<0) ad=0;
            //else ad = 180;
            //theta = ad - 180.0/M_PI*Foam::atan(y/x);
            theta = ad + pm*180.0/M_PI*Foam::atan((y-2*r)/x);
            std::cout<<100*Cf[facei].component(0) << "\t" << 100*Cf[facei].component(1)-2.0;
            Info<<"\t"<< theta;
            Info<<"\t"<<(q[facei]&Sf[facei])/mag(Sf[facei]);
            Info<<"\t"<<(F[facei]&Sf[facei])/mag(Sf[facei]);
            Info<<"\t"<<magSqr(F[facei]- (F[facei]&Sf[facei])*Sf[facei]/magSqr(Sf[facei])) << endl;
        }
    }
    return 0;
}

// ************************************************************************* //
