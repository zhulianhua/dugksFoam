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

\*---------------------------------------------------------------------------*/

#include "fvDVM.H"
#include "discreteVelocity.H"
#include "constants.H"
#include "fixedGradientFvPatchField.H"
#include "calculatedMaxwellFvPatchField.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;


//const Foam::word
//Foam::radiation::radiativeIntensityRay::intensityPrefix("ILambda");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::discreteVelocity::discreteVelocity
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
)
:
    dvm_(dvm),
    mesh_(mesh),
    time_(time),
    weight_(weight),
    xi_(xi),
    myDVid_(DVid),
    symXtargetDVid_(symXtargetDVid),
    symYtargetDVid_(symYtargetDVid),
    symZtargetDVid_(symZtargetDVid),
    gTildeVol_
    (
        IOobject
        (
            "gTildeVol" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass*pow3(dimTime/dimLength)/pow3(dimLength), 0.0
        )
    ),
    hTildeVol_
    (
        IOobject
        (
            "hTildeVol" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass*(dimTime/dimLength)/pow3(dimLength), 0.0
        )
    ),
    gBarPvol_
    (
        IOobject
        (
            "gBarPvol" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass*pow3(dimTime/dimLength)/pow3(dimLength), 0.0
        ),
        "fixedGradient"
    ),
    hBarPvol_
    (
        IOobject
        (
            "hBarPvol" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass*(dimTime/dimLength)/pow3(dimLength), 0.0
        ),
        "fixedGradient"
    ),
    gSurf_
    (
        IOobject
        (
            "gSurf" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass*pow3(dimTime/dimLength)/pow3(dimLength), 0.0
        )
    ),
    hSurf_
    (
        IOobject
        (
            "hSurf" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "0", dimMass*(dimTime/dimLength)/pow3(dimLength), 0.0
        )
    ),
    gBarPgrad_
    (
        IOobject
        (
            "gBarPgrad" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "0", gBarPvol_.dimensions()/dimLength, vector(0,0,0)
        ),
        "zeroGradient"
    ),
    hBarPgrad_
    (
        IOobject
        (
            "hBarPgrad" + name(DVid),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "0", hBarPvol_.dimensions()/dimLength, vector(0,0,0)
        ),
        "zeroGradient"
    )
{
    initDFtoEq();
    setBCtype();
    initBoundaryField();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::discreteVelocity::~discreteVelocity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::discreteVelocity::initDFtoEq()
{
    GeometricField<vector, fvPatchField, volMesh> qVolIni
    (
        IOobject
        (
            "qVolIni",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "0", dvm_.qVol().dimensions(), vector(0,0,0)
        )
    );

    equilibriumShakhov
    (
        gTildeVol_,
        hTildeVol_,
        dvm_.rhoVol(), 
        dvm_.Uvol(), 
        dvm_.Tvol(), 
        qVolIni
    );

}

void Foam::discreteVelocity::setBCtype()
{
    // Only from rho's B.C can we determine surfaceScalarField f's B.C.
    // for all patchi fo h/g barP, set to zeroGradient
    // May improve for inlet ingoing DF

    const GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = dvm_.rhoVol().boundaryField();

    forAll(rhoBCs, patchi)
    {
        if (rhoBCs[patchi].type() == "fixedValue") //inlet
        {
            gSurf_.boundaryField().set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    "mixed", mesh_.boundary()[patchi], gSurf_
                )
            );
            hSurf_.boundaryField().set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    "mixed", mesh_.boundary()[patchi], hSurf_
                )
            );
        }
        else if (rhoBCs[patchi].type() == "zeroGradient") //outlet
        {
            gSurf_.boundaryField().set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    "zeroGradient", mesh_.boundary()[patchi], gSurf_
                )
            );
            hSurf_.boundaryField().set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    "zeroGradient", mesh_.boundary()[patchi], hSurf_
                )
            );
        }
        else if (rhoBCs[patchi].type() == "calculatedMaxwell") //maxwellWall
        {
            gSurf_.boundaryField().set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    "maxwellWall", mesh_.boundary()[patchi], gSurf_
                )
            );
            hSurf_.boundaryField().set
            (
                patchi, 
                fvsPatchField<scalar>::New
                (
                    "maxwellWall", mesh_.boundary()[patchi], hSurf_
                )
            );
        }
    }
}

void Foam::discreteVelocity::initBoundaryField()
{
    // fvsPatchField gSurf_ hSurf_  with "mixed type" init field
    // Here we only consider g but not h, since the eq of h is zero, which is 
    // the default value of mixedFvsPatchField
    // Bug: h is not zero
    GeometricField<scalar, fvsPatchField, surfaceMesh>::GeometricBoundaryField& 
        gBCs = gSurf_.boundaryField();
    GeometricField<scalar, fvsPatchField, surfaceMesh>::GeometricBoundaryField& 
        hBCs = hSurf_.boundaryField();


    forAll(gBCs, patchi)
    {
        if(gBCs[patchi].type() == "mixed" )
        {
            equilibriumMaxwell
            (
                gBCs[patchi],
                hBCs[patchi],
                dvm_.rhoVol().boundaryField()[patchi],
                dvm_.Uvol().boundaryField()[patchi],
                dvm_.Tvol().boundaryField()[patchi]
            );
        }
    }
}

void Foam::discreteVelocity::updateGHbarPvol()
{
    //- get delta t
    dimensionedScalar dt = time_.deltaT();
    volScalarField gEq
    (
        IOobject
        (
            "gEq",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", gTildeVol_.dimensions(), 0)
    );

    volScalarField hEq
    (
        IOobject
        (
            "hEq",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", hTildeVol_.dimensions(), 0)
    );

    volScalarField relaxFactor
    (
        IOobject
        (
            "relaxFactor",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0) // dimLess
    );

    //- get relaxtion factor 
    relaxFactor = 1.5*dt/(2.0*dvm_.tauVol() + dt);

    //- get gEq and hEq
    equilibriumShakhov
    (
        gEq,
        hEq, 
        dvm_.rhoVol(), 
        dvm_.Uvol(), 
        dvm_.Tvol(), 
        dvm_.qVol() 
    );
    gBarPvol_ = (1.0 - relaxFactor)*gTildeVol_ + relaxFactor*gEq;
    hBarPvol_ = (1.0 - relaxFactor)*hTildeVol_ + relaxFactor*hEq;

    gBarPvol_.correctBoundaryConditions(); // NOTE: check if the newly defined zeroGradientFvsPatchField 
    hBarPvol_.correctBoundaryConditions();
}

void Foam::discreteVelocity::updateGHbarSurf()
{

  // Pre-setup
    // 1. correct the boundary value of gBarPvol_
    //    the gradient at boundary is known
    // 2. get the gradient
    //
    gBarPgrad_ = fvc::grad(gBarPvol_); 
    hBarPgrad_ = fvc::grad(hBarPvol_);
    // The DVMsymmetry is  rocessed automatically in fvc::grad operator

    //
    // 3. correct the boundary value of the grad field
    //    to be used at next time
    gBarPgrad_.correctBoundaryConditions();
    hBarPgrad_.correctBoundaryConditions();

    // 4. patch the normal component of boundary value of the grad 
    //    to the gradient field of the fixed 
    //    gradient feild of the gBarPvol_ and 
    //    hBarPvol_ ...
    //    NOTE: we need the surfaceNormal Gradient

    forAll(gBarPgrad_.boundaryField(), patchi)
    {
        const vectorField n
        (
            mesh_.Sf().boundaryField()[patchi]
           /mesh_.magSf().boundaryField()[patchi]
        );

        if ( 
               gBarPvol_.boundaryField()[patchi].type() != "empty" 
            && gBarPvol_.boundaryField()[patchi].type() != "processor"
            && gBarPvol_.boundaryField()[patchi].type() != "symmetryPlane"
           ) // only for fixed gradient g/hBarPvol
        {
            // normal component of the grad field
            fixedGradientFvPatchField<scalar>& gBarPvolPatch = 
                refCast<fixedGradientFvPatchField<scalar> >
                (gBarPvol_.boundaryField()[patchi]);

            fixedGradientFvPatchField<scalar>& hBarPvolPatch = 
                refCast<fixedGradientFvPatchField<scalar> >
                (hBarPvol_.boundaryField()[patchi]);

            forAll(gBarPvolPatch, pFacei)
            {
                gBarPvolPatch.gradient()[pFacei] =
                    gBarPgrad_.boundaryField()[patchi][pFacei]&n[pFacei];
                hBarPvolPatch.gradient()[pFacei] =
                    hBarPgrad_.boundaryField()[patchi][pFacei]&n[pFacei];
            }
        }
    }

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField& Cf = mesh_.Cf();
    const vectorField& Sf = mesh_.Sf();
    const vectorField& C = mesh_.C();

    const Field<scalar>& iGbarPvol = gBarPvol_;
    const Field<scalar>& iHbarPvol = hBarPvol_;
    const Field<vector>& iGbarPgrad = gBarPgrad_;
    const Field<vector>& iHbarPgrad = hBarPgrad_;

    // This is what we want to update in this function
    Field<scalar>& iGsurf = gSurf_;
    Field<scalar>& iHsurf = hSurf_;

    scalar dt = time_.deltaTValue();
    vector xii = xi_.value();

    // internal faces first
    forAll(owner, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];
        if ((xii&Sf[facei]) >=  0) // comming from own
        {

            iGsurf[facei] = iGbarPvol[own] 
              + (iGbarPgrad[own]&(Cf[facei] - C[own] - 0.5*xii*dt));

            iHsurf[facei] = iHbarPvol[own]
              + (iHbarPgrad[own]&(Cf[facei] - C[own] - 0.5*xii*dt));

        }
        // Debug, no = 0, =0 put to > 0
        else// if ((xii&Sf[facei]) < 0) // comming form nei
        {
            iGsurf[facei] = iGbarPvol[nei]
              + (iGbarPgrad[nei]&(Cf[facei] - C[nei] - 0.5*xii*dt));
            iHsurf[facei] = iHbarPvol[nei]
              + (iHbarPgrad[nei]&(Cf[facei] - C[nei] - 0.5*xii*dt));
        }
        /*
         *else // NOTE: may be improved by using weighted combination
         *{
         *    iGsurf[facei] = 0.5*
         *    (
         *        iGbarPvol[nei] + ((iGbarPgrad[nei])
         *       &(Cf[facei] - C[nei] - xii*dt))
         *      + iGbarPvol[own] + ((iGbarPgrad[own])
         *       &(Cf[facei] - C[nei] - xii*dt))
         *    );
         *    iHsurf[facei] = 0.5*
         *    (
         *        iHbarPvol[own] + ((iHbarPgrad[own])
         *       &(Cf[facei] - C[nei] - xii*dt))
         *      + iHbarPvol[nei] + ((iHbarPgrad[nei])
         *       &(Cf[facei] - C[nei] - xii*dt))
         *    );
         *}
         */
    }

  // boundary faces
    forAll(gSurf_.boundaryField(), patchi)
    {
        word type = gSurf_.boundaryField()[patchi].type();
        fvsPatchField<scalar>& gSurfPatch = gSurf_.boundaryField()[patchi];
        fvsPatchField<scalar>& hSurfPatch = hSurf_.boundaryField()[patchi];
        const fvsPatchField<vector>& SfPatch =
            mesh_.Sf().boundaryField()[patchi];
        const fvsPatchField<vector>& CfPatch =
            mesh_.Cf().boundaryField()[patchi];
        const labelUList& faceCells = mesh_.boundary()[patchi].faceCells();
        
        //- NOTE: outging DF can be treate unifily for all BCs, including processor BC
        if (type == "zeroGradient")
        {
            gSurfPatch == gBarPvol_.boundaryField()[patchi].patchInternalField();
            hSurfPatch == hBarPvol_.boundaryField()[patchi].patchInternalField();
        }
        else if (type == "mixed")
        {
            //check each boundary face in the patch
            forAll(gSurfPatch, facei)
            {
                //out or in ?
                if ((xii&SfPatch[facei]) > 0 ) // outgoing
                {
                    gSurfPatch[facei] = iGbarPvol[faceCells[facei]] 
                      + ((iGbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - xii*dt));
                    hSurfPatch[facei] = iHbarPvol[faceCells[facei]] 
                      + ((iHbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - xii*dt));
                //incoming and parallel to face, not changed.
                }
            }
        }
        else if (type == "maxwellWall")
        {
            calculatedMaxwellFvPatchField<scalar>& rhoPatch = 
                refCast<calculatedMaxwellFvPatchField<scalar> >
                (dvm_.rhoVol().boundaryField()[patchi]); //DEBUG

            forAll(gSurfPatch, facei)
            {
                vector faceSf= SfPatch[facei];
                if ((xii&faceSf) >  0 ) // outgoing
                {
                    gSurfPatch[facei] = iGbarPvol[faceCells[facei]] 
                      + ((iGbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - xii*dt));
                    hSurfPatch[facei] = iHbarPvol[faceCells[facei]] 
                      + ((iHbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - xii*dt));

                    rhoPatch.outGoing()[facei] += //add outgoing normal momentum flux to outGoing container
                        weight_*(xii&faceSf)*gSurfPatch[facei];
                }
            }
        }
        else if (type == "processor") // parallel
        {
            forAll(gSurfPatch, facei)
            {
                vector faceSf= SfPatch[facei];
                if ((xii&faceSf) >  0 ) // outgoing
                {
                    gSurfPatch[facei] = iGbarPvol[faceCells[facei]] 
                      + ((iGbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - xii*dt));
                    hSurfPatch[facei] = iHbarPvol[faceCells[facei]] 
                      + ((iHbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - xii*dt));

                } 
                else //incomming from processor boundaryField
                {
                    gSurfPatch[facei] = gBarPvol_.boundaryField()[patchi][facei]
                      + ((gBarPgrad_.boundaryField()[patchi][facei])
                       &(CfPatch[facei] - mesh_.C().boundaryField()[patchi][facei] - xii*dt));
                    hSurfPatch[facei] = hBarPvol_.boundaryField()[patchi][facei]
                      + ((hBarPgrad_.boundaryField()[patchi][facei])
                       &(CfPatch[facei] - mesh_.C().boundaryField()[patchi][facei] - xii*dt));
                }

            }
        }
        else if (type == "symmetryPlane")
        {
            forAll(gSurfPatch, facei)
            {
                vector faceSf= SfPatch[facei];
                if ((xii&faceSf) >  0 ) // outgoing, incomming shoud be proceed after this function
                {
                    gSurfPatch[facei] = iGbarPvol[faceCells[facei]] 
                      + ((iGbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - xii*dt));
                    hSurfPatch[facei] = iHbarPvol[faceCells[facei]] 
                      + ((iHbarPgrad[faceCells[facei]])
                       &(CfPatch[facei] - C[faceCells[facei]] - xii*dt));
                } 
            }
        }
    }
}

void Foam::discreteVelocity::updateGHbarSurfMaxwellWallIn()
{
    vector xii = xi_.value();
    forAll(gSurf_.boundaryField(), patchi)
    {
        if (gSurf_.boundaryField()[patchi].type() == "maxwellWall")
        {
            fvsPatchScalarField& gSurfPatch = gSurf_.boundaryField()[patchi];
            fvsPatchScalarField& hSurfPatch = hSurf_.boundaryField()[patchi];
            const fvPatchScalarField& rhoVolPatch = 
                dvm_.rhoVol().boundaryField()[patchi];
            const fvPatchVectorField& UvolPatch = 
                dvm_.Uvol().boundaryField()[patchi];
            const fvPatchScalarField& TvolPatch = 
                dvm_.Tvol().boundaryField()[patchi];
            const fvsPatchVectorField& SfPatch = 
                mesh_.Sf().boundaryField()[patchi];
            forAll(gSurfPatch, facei)
            {
                vector faceSf = SfPatch[facei];
                if ((xii & faceSf) <= 0) // incomming
                {
                    // set to maxwellian
                    gSurfPatch[facei] = rhoVolPatch[facei]
                       *equilibriumMaxwellByRho
                        (
                            UvolPatch[facei],
                            TvolPatch[facei]
                        );

                    //set hSurf at maxwellWall to zero! , WRONG!!!
                    hSurfPatch[facei] = 
                        gSurfPatch[facei]*(dvm_.R().value()*TvolPatch[facei])
                       *(dvm_.KInner() + 3 - mesh_.nSolutionD());
                }
            }
        }
    }
}

void Foam::discreteVelocity::updateGHbarSurfSymmetryIn()
{
    forAll(gSurf_.boundaryField(), patchi)
    {
        if (gSurf_.boundaryField()[patchi].type() == "symmetryPlane")
        {
            fvsPatchScalarField& gSurfPatch = gSurf_.boundaryField()[patchi];
            fvsPatchScalarField& hSurfPatch = hSurf_.boundaryField()[patchi];
            vector faceSf = mesh_.Sf().boundaryField()[patchi][0];
            vector nomlizedDirec = faceSf/mag(faceSf);
            label targetDVid = nomlizedDirec
                & vector(symXtargetDVid_, symYtargetDVid_, symZtargetDVid_);
            forAll(gSurfPatch, facei)
            {
                gSurfPatch[facei] = dvm_.DVi(targetDVid).gSurf().boundaryField()
                    [patchi][facei];
                hSurfPatch[facei] = dvm_.DVi(targetDVid).hSurf().boundaryField()
                    [patchi][facei];
            }
        }
    }
}

void Foam::discreteVelocity::updateGHsurf()
{
    //- get delta t
    dimensionedScalar h = 0.5*time_.deltaT();

    surfaceScalarField gEq
    (
        IOobject
        (
            "gEq",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", gTildeVol_.dimensions(), 0)
    );

    surfaceScalarField hEq
    (
        IOobject
        (
            "hEq",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", hTildeVol_.dimensions(), 0)
    );

    surfaceScalarField relaxFactor
    (
        IOobject
        (
            "relaxFactor",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        scalar(0)
    );

    //- get relaxtion factor 
    relaxFactor = h/(2*dvm_.tauSurf() + h);

    //- get gEq and hEq
    equilibriumShakhov
    (
        gEq,
        hEq, 
        dvm_.rhoSurf(), 
        dvm_.Usurf(), 
        dvm_.Tsurf(), 
        dvm_.qSurf() 
    );

    gSurf_ = (1.0 - relaxFactor)*gSurf_ + relaxFactor*gEq;
    hSurf_ = (1.0 - relaxFactor)*hSurf_ + relaxFactor*hEq;
    // NOTE: here the boundar face value are not computed

    vector xii = xi_.value();
    // We do it mannuly for the outgoing DF?? Why??
    forAll(gSurf_.boundaryField(), patchi)
    {
        fvsPatchScalarField& gSurfPatch = gSurf_.boundaryField()[patchi];
        fvsPatchScalarField& hSurfPatch = hSurf_.boundaryField()[patchi];
        fvsPatchScalarField& gEqPatch   =    gEq.boundaryField()[patchi];
        fvsPatchScalarField& hEqPatch   =    hEq.boundaryField()[patchi];

        const fvsPatchVectorField& SfPatch = mesh_.Sf().boundaryField()[patchi];

        fvsPatchScalarField& relaxFactorPatch = 
            relaxFactor.boundaryField()[patchi];

        forAll(gSurfPatch, facei)
        {
            if ( (xii&(SfPatch[facei])) > 0  ) // Here, if delted , the free stream BC may bo s problem
            {
                gSurfPatch[facei] = (1.0 - relaxFactorPatch[facei])
                    *gSurfPatch[facei]
                   + relaxFactorPatch[facei]*gEqPatch[facei];
                hSurfPatch[facei] = (1.0 - relaxFactorPatch[facei])
                    *hSurfPatch[facei]
                   + relaxFactorPatch[facei]*hEqPatch[facei];
            }
        }
    }
}

void Foam::discreteVelocity::updateGHtildeVol()
{
    // store the gTildePlus in gTilde
    gTildeVol_ = -1.0/3*gTildeVol_ + 4.0/3*gBarPvol_;
    hTildeVol_ = -1.0/3*hTildeVol_ + 4.0/3*hBarPvol_;

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();
    const vectorField Sf = mesh_.Sf();
    const scalarField V = mesh_.V();
    const scalar dt = time_.deltaTValue();
    const vector xii = xi_.value();

    // internal faces
    forAll(owner, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];
        gTildeVol_[own] -= ((xii&Sf[facei])*gSurf_[facei]*dt/V[own]);
        gTildeVol_[nei] += ((xii&Sf[facei])*gSurf_[facei]*dt/V[nei]);
        hTildeVol_[own] -= ((xii&Sf[facei])*hSurf_[facei]*dt/V[own]);
        hTildeVol_[nei] += ((xii&Sf[facei])*hSurf_[facei]*dt/V[nei]);
    }

    // boundary faces
    forAll(gSurf_.boundaryField(), patchi)
    {
        const fvsPatchField<scalar>& gSurfPatch =
            gSurf_.boundaryField()[patchi];
        const fvsPatchField<scalar>& hSurfPatch =
            hSurf_.boundaryField()[patchi];
        const fvsPatchField<vector>& SfPatch =
            mesh_.Sf().boundaryField()[patchi];
        const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
        forAll(pOwner, pFacei)
        {
            const label own = pOwner[pFacei];
            gTildeVol_[own] -= (xii&SfPatch[pFacei]) 
               *gSurfPatch[pFacei]*dt/V[own];
            hTildeVol_[own] -= (xii&SfPatch[pFacei]) 
               *hSurfPatch[pFacei]*dt/V[own];
        }
    }

}

template <template<class> class PatchType, class GeoMesh> 
Foam::tmp<Foam::GeometricField<Foam::scalar, PatchType, GeoMesh> >
Foam::discreteVelocity::equilibriumMaxwell
(
    const GeometricField<scalar, PatchType, GeoMesh>& rho,
    const GeometricField<vector, PatchType, GeoMesh>& U,
    const GeometricField<scalar, PatchType, GeoMesh>& T
)
{
    const dimensionedScalar R = dvm_.R();
    const label D = mesh_.nSolutionD();
    const dimensionedScalar vUnit("vUnit", dimLength/dimTime, 1);

    tmp<GeometricField<scalar, PatchType, GeoMesh> > tEqu
    (
        new GeometricField<scalar, PatchType, GeoMesh>
        (
            IOobject
            (
                "equ",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar( "0", gTildeVol_.dimensions(), 0)
        )
    );
    GeometricField<scalar, PatchType, GeoMesh>& equ = tEqu();
    equ = rho/pow(sqrt(2.0*pi*R*T),D)*exp(-magSqr(U - xi_)/(2.0*R*T))
        /pow(vUnit, 3-D);
    return tEqu;
}

template <template<class> class PatchType, class GeoMesh> 
void Foam::discreteVelocity::equilibriumShakhov
(
    GeometricField<scalar, PatchType, GeoMesh>& gEq,
    GeometricField<scalar, PatchType, GeoMesh>& hEq,
    const GeometricField<scalar, PatchType, GeoMesh>& rho,
    const GeometricField<vector, PatchType, GeoMesh>& U,
    const GeometricField<scalar, PatchType, GeoMesh>& T,
    const GeometricField<vector, PatchType, GeoMesh>& q
)
{
    dimensionedScalar R = dvm_.R();
    dimensionedScalar Pr = dvm_.Pr();
    label D = mesh_.nSolutionD();
    label K = dvm_.KInner();
    dimensionedScalar vUnit("vUnit", dimLength/dimTime, 1);

    GeometricField<scalar, PatchType, GeoMesh> cSqrByRT 
        = magSqr(U - xi_)/(R*T);

    GeometricField<scalar, PatchType, GeoMesh> cqBy5pRT 
        = ((xi_ - U)&q)/(5.0*rho*R*T*R*T);

    GeometricField<scalar, PatchType, GeoMesh> gEqBGK 
        = rho/pow(sqrt(2.0*pi*R*T),D)*exp(-cSqrByRT/2.0)/pow(vUnit, 3-D);

    gEq = ( 1.0 + (1.0 - Pr)*cqBy5pRT*(cSqrByRT - D - 2.0) )*gEqBGK;
    hEq = ( (K + 3.0 - D) + (1.0 - Pr)*cqBy5pRT*((cSqrByRT - D)*(K + 3.0 - D) - 2*K) )*gEqBGK*R*T;
}

void Foam::discreteVelocity::equilibriumMaxwell
(
    Foam::fvsPatchScalarField& geq,
    Foam::fvsPatchScalarField& heq,
    const Foam::fvPatchScalarField&  rho,
    const Foam::fvPatchVectorField&    U,
    const Foam::fvPatchScalarField&    T
)
{
    //No need to consider the dimensionSet!!
    scalar Ri = dvm_.R().value(); //but R has dimensionSet
    //vector xii = xi_;
    label D = mesh_.nSolutionD();
    vector xii = xi_.value();
    geq == rho/pow(sqrt(2.0*pi*Ri*T),D)*exp(-magSqr(U - xii)/(2.0*Ri*T));
    heq == (dvm_.KInner() + 3-D)*Ri*T*geq;
}


Foam::scalar Foam::discreteVelocity::equilibriumMaxwellByRho
(
 const  vector U,
 const  scalar T
)
{
    scalar feqByRho;
    scalar Ri = dvm_.R().value(); //but R has dimensionSet
    label D = mesh_.nSolutionD();
    vector xii = xi_.value();
    feqByRho = 1.0/pow(sqrt(2.0*pi*Ri*T),D)*exp(-magSqr(U - xii)/(2.0*Ri*T));
    return feqByRho;
}

// ************************************************************************* //
