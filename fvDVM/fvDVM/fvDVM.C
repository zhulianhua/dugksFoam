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

#include "fvDVM.H"
#include "constants.H"
#include "fvm.H"
#include "calculatedMaxwellFvPatchField.H"
#include "scalarIOList.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvDVM, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvDVM::setDVgrid
(
    scalarField& weights,
    scalarField& Xis,
    scalar xiMin,
    scalar xiMax,
    label nXi
)
{
    // Read from file ./constant/Xis and ./constant/weights
    Info << "Reading xi list" << endl;
    scalarIOList xiList
    (
        IOobject
        (
             "Xis",
             time_.time().constant(),
             mesh_,
             IOobject::MUST_READ,
             IOobject::NO_WRITE
        )
    );

    Info << "Reading weight list" << endl;
    scalarIOList weightList
    (
        IOobject
        (
             "weights",
             time_.time().constant(),
             mesh_,
             IOobject::MUST_READ,
             IOobject::NO_WRITE
        )
    );

    // Check
    if (
            xiList.size() != nXi 
         || weightList.size() != nXi
       )
    {
        WarningIn("setDVgrid")
            << "Num of discreteVelocity not consistent" << endl;
        std::exit(-1);
    }


    if ( 
           xiList[0]  != xiMin
        || xiList[nXi-1] != xiMax)
    {
        WarningIn("setDVgrid")
            << "xi not consistant" <<endl;
        std::exit(-1);
    }


    for (label i = 0; i < nXi ; i++)
    {
        weights[i] = weightList[i];
        Xis[i] = xiList[i];
    }
}

void Foam::fvDVM::initialiseDV()
{
    scalarField weights1D(nXiPerDim_);
    scalarField Xis(nXiPerDim_);

    //get discrete velocity points and weights
    setDVgrid
    (
         weights1D,
         Xis, 
         xiMin_.value(), 
         xiMax_.value(), 
         nXiPerDim_
    );

    if (mesh_.nSolutionD() == 3)    //3D(X & Y & Z)
    {
        nXiX_ = nXiY_ = nXiZ_ = nXiPerDim_;
        nXi_ = nXiX_*nXiY_*nXiZ_;
        DV_.setSize(nXi_);
        label i = 0;
        for (label iz = 0; iz < nXiZ_; iz++)
        {
            for (label iy = 0; iy < nXiY_; iy++)
            {
                for (label ix = 0; ix < nXiZ_; ix++)
                {
                    scalar weight = weights1D[iz]*weights1D[iy]*weights1D[ix];
                    dimensionedVector xi
                    (
                        "xi", 
                        dimLength/dimTime, 
                        vector(Xis[ix], Xis[iy], Xis[iz])
                    );

                    label symXtargetDVid = iz*nXiY_*nXiX_ + iy*nXiX_ + (nXiX_ - ix -1);
                    label symYtargetDVid = iz*nXiY_*nXiX_ + (nXiY_ - iy - 1)*nXiX_ + ix;
                    label symZtargetDVid = (nXiZ_ - iz -1)*nXiY_*nXiX_ + iy*nXiX_ + ix;
                    DV_.set
                    (
                        i,
                        new discreteVelocity
                        (
                            *this,
                            mesh_,
                            time_,
                            weight,
                            xi,
                            i,
                            symXtargetDVid,
                            symYtargetDVid,
                            symZtargetDVid
                        )
                    );
                    i++;
                }
            }
        }
    }
    else
    {
        if (mesh_.nSolutionD() == 2)    //2D (X & Y)
        {
            nXiX_ = nXiY_ = nXiPerDim_;
            nXiZ_ = 1;
            nXi_ = nXiX_*nXiY_*nXiZ_;
            DV_.setSize(nXi_);
            label i = 0;
            for (label iy = 0; iy < nXiY_; iy++)
            {
                for (label ix = 0; ix < nXiX_; ix++)
                {
                    scalar weight = weights1D[iy]*weights1D[ix]*1;
                    dimensionedVector xi
                    (
                        "xi",
                        dimLength/dimTime, 
                        vector(Xis[ix], Xis[iy], 0.0)
                    );
                    label symXtargetDVid = iy*nXiX_ + (nXiX_ - ix -1);
                    label symYtargetDVid = (nXiY_ - iy - 1)*nXiX_ + ix;
                    label symZtargetDVid = 0;
                    DV_.set
                    (
                        i,
                        new discreteVelocity
                        (
                            *this,
                            mesh_,
                            time_,
                            weight,
                            xi,
                            i,
                            symXtargetDVid,
                            symYtargetDVid,
                            symZtargetDVid
                        )
                    );
                    i++;
                }
            }
        }
        else    //1D (X)
        {
            nXiX_ = nXiPerDim_;
            nXiY_ = nXiZ_ = 1;
            nXi_ = nXiX_*nXiY_*nXiZ_;
            DV_.setSize(nXi_);
            label i = 0;
            for (label ix = 0; ix < nXiX_; ix++)
            {
                scalar weight = weights1D[ix]*1*1;
                dimensionedVector xi
                (
                    "xi",
                    dimLength/dimTime,
                    vector(Xis[ix], 0.0, 0.0)
                );
                label symXtargetDVid = (nXiX_ - ix -1);
                label symYtargetDVid = 0;
                label symZtargetDVid = 0;
                DV_.set
                (
                    i,
                    new discreteVelocity
                    (
                        *this,
                        mesh_,
                        time_,
                        weight,
                        xi,
                        i,
                        symXtargetDVid,
                        symYtargetDVid,
                        symZtargetDVid
                    )
                );
                i++;
            }
        }
    }

    Info<< "fvDVM : Allocated " << DV_.size()
        << " discrete velocities" << endl;
}


void Foam::fvDVM::setCalculatedMaxwellRhoBC()
{
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
    forAll(rhoBCs, patchi)
    {
        if (rhoBCs[patchi].type() == "calculatedMaxwell")
        {

            const vectorField& SfPatch = mesh_.Sf().boundaryField()[patchi];
            calculatedMaxwellFvPatchField<scalar>& rhoPatch = 
                refCast<calculatedMaxwellFvPatchField<scalar> >(rhoBCs[patchi]);
            const vectorField& Upatch = Uvol_.boundaryField()[patchi];
            const scalarField& Tpatch = Tvol_.boundaryField()[patchi];
            forAll(rhoPatch, facei)
            {
                vector faceSf = SfPatch[facei];
                rhoPatch.inComingByRho()[facei] = 0; // set to zero
                forAll(DV_, dvi) // add one by one
                {
                    vector xi = DV_[dvi].xi().value();
                    scalar weight = DV_[dvi].weight();
                    if ( (xi & faceSf) < 0) //inComing
                    {
                        rhoPatch.inComingByRho()[facei] += 
                          - weight*(xi & faceSf)
                          *DV_[dvi].equilibriumMaxwellByRho
                          (
                              Upatch[facei], 
                              Tpatch[facei]
                          );
                    }
                }
            }
        }
    }
}


void Foam::fvDVM::updateGHbarPvol()
{
    forAll(DV_, DVid)
        DV_[DVid].updateGHbarPvol();
}


void Foam::fvDVM::updateGHbarSurf()
{
    forAll(DV_, DVid)
        DV_[DVid].updateGHbarSurf();
}


void Foam::fvDVM::updateMaxwellWallRho()
{
    rhoVol_.correctBoundaryConditions(); //Bug here !!!
}


void Foam::fvDVM::updateGHbarSurfMaxwellWallIn()
{
    forAll(DV_, DVid)
        DV_[DVid].updateGHbarSurfMaxwellWallIn();
}


void Foam::fvDVM::updateMacroSurf()
{
    // Init to zero before add one DV by one DV
    rhoSurf_ =  dimensionedScalar("0", rhoSurf_.dimensions(), 0);
    Usurf_ = dimensionedVector("0", Usurf_.dimensions(), vector(0, 0, 0));
    Tsurf_ = dimensionedScalar("0", Tsurf_.dimensions(), 0);
    qSurf_ = dimensionedVector("0", qSurf_.dimensions(), vector(0, 0, 0));
    stressSurf_ = dimensionedTensor
        (
            "0", 
            stressSurf_.dimensions(), 
            pTraits<tensor>::zero
        );

    // Conserved variable, now zero as the prime variable 
    // has been set to zero
    surfaceVectorField rhoUsurf = rhoSurf_*Usurf_;
    surfaceScalarField rhoEsurf = rhoSurf_*magSqr(Usurf_);

    forAll(DV_, dvi)
    {
        discreteVelocity& dv = DV_[dvi];
        rhoSurf_ += dXiCellSize_*dv.weight()*dv.gSurf();
        rhoUsurf += dXiCellSize_*dv.weight()*dv.gSurf()*dv.xi();
        rhoEsurf += 0.5*dXiCellSize_*dv.weight()
           *(
                dv.gSurf()*magSqr(dv.xi()) 
              + dv.hSurf()
            );
    }

    //- get Prim. from Consv.

    Usurf_ = rhoUsurf/rhoSurf_;
    Tsurf_ = (rhoEsurf - 0.5*rhoSurf_*magSqr(Usurf_))/((KInner_ + 3)/2.0*R_*rhoSurf_);

    //- update tau
    tauSurf_ = updateTau(Tsurf_, rhoSurf_);

    //- peculiar vel.
    surfaceVectorField c = Usurf_;

    //-get heat flux
    forAll(DV_, dvi)
    {
        discreteVelocity& dv = DV_[dvi];
        c = dv.xi() - Usurf_;
        qSurf_ += 0.5*dXiCellSize_*dv.weight()*c
            *(
                 magSqr(c)*dv.gSurf() 
               + dv.hSurf()
             );
        stressSurf_ += 
            dXiCellSize_*dv.weight()*dv.gSurf()*c*c;
    }

    //- correction for bar to original
    qSurf_ = 2.0*tauSurf_/(2.0*tauSurf_ + 0.5*time_.deltaT()*Pr_)*qSurf_;
    stressSurf_ = 
        2.0*tauSurf_/(2.0*tauSurf_ + 0.5*time_.deltaT())*stressSurf_;
}


void Foam::fvDVM::updateGHsurf()
{
    forAll(DV_, DVid)
        DV_[DVid].updateGHsurf();
}


void Foam::fvDVM::updateGHtildeVol()
{
    forAll(DV_, DVid)
        DV_[DVid].updateGHtildeVol();
}

void Foam::fvDVM::updateMacroVol()
{

    volVectorField rhoUvol = rhoVol_*Uvol_;
    volScalarField rhoEvol = rhoVol_*(magSqr(Uvol_) + (KInner_ + 3)/2.0*R_*Tvol_);

    qVol_ = dimensionedVector("0", qVol_.dimensions(), vector(0, 0, 0));

    if(macroFlux_ == "no") // update cell macro by moment from DF
    {
        rhoVol_ =  dimensionedScalar("0", rhoVol_.dimensions(), 0);
        rhoUvol = dimensionedVector("0", rhoUvol.dimensions(), vector(0, 0, 0));
        rhoEvol = dimensionedScalar("0", rhoEvol.dimensions(), 0);

        forAll(DV_, dvi)
        {
            discreteVelocity& dv = DV_[dvi];
            rhoVol_ += dXiCellSize_*dv.weight()*dv.gTildeVol();
            rhoUvol += dXiCellSize_*dv.weight()*dv.gTildeVol()*dv.xi();
            rhoEvol += 0.5*dXiCellSize_*dv.weight()
               *(
                    magSqr(dv.xi())*dv.gTildeVol() 
                  + dv.hTildeVol()
                );
        }
    }
    else // update by macro flux
    {
        const labelUList& owner = mesh_.owner();
        const labelUList& neighbour = mesh_.neighbour();
        const vectorField Sf = mesh_.Sf();
        const scalarField V = mesh_.V();
        const scalar dt = time_.deltaTValue();

        // internal faces
        forAll(owner, facei)
        {
            const label own = owner[facei];
            const label nei = neighbour[facei];
            rhoVol_[own] -= (rhoSurf_[facei]*Usurf_[facei]&Sf[facei])*dt/V[own];
            rhoVol_[nei] += (rhoSurf_[facei]*Usurf_[facei]&Sf[facei])*dt/V[nei];
            rhoUvol[own] -= (rhoSurf_[facei]*Usurf_[facei]*Usurf_[facei]
                            + stressSurf_[facei])&Sf[facei]*dt/V[own];
            rhoUvol[nei] += (rhoSurf_[facei]*Usurf_[facei]*Usurf_[facei]
                            + stressSurf_[facei])&Sf[facei]*dt/V[nei];
            scalar rhoEsurf = 
                rhoSurf_[facei]
               *(magSqr(Usurf_[facei]) + (KInner_ + 3)/2.0*R_.value()*Tsurf_[facei]);

            rhoEvol[own] -= (rhoEsurf*Usurf_[facei] + qSurf_[facei])
                            &Sf[facei]*dt/V[own];
            rhoEvol[nei] += (rhoEsurf*Usurf_[facei] + qSurf_[facei])
                            &Sf[facei]*dt/V[nei];
        }
        // boundary faces
        forAll(rhoSurf_.boundaryField(), patchi)
        {
            const fvsPatchField<scalar>& rhoSurfPatch =
                rhoSurf_.boundaryField()[patchi];
            const fvsPatchField<vector>& UsurfPatch =
                Usurf_.boundaryField()[patchi];
            const fvsPatchField<scalar>& TsurfPatch =
                Tsurf_.boundaryField()[patchi];
            const fvsPatchField<tensor>& stressSurfPatch =
                stressSurf_.boundaryField()[patchi];
            const fvsPatchField<vector>& qSurfPatch =
                qSurf_.boundaryField()[patchi];
            const fvsPatchField<vector>& SfPatch =
                mesh_.Sf().boundaryField()[patchi];

            const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
            forAll(pOwner, pFacei)
            {
                const label own = pOwner[pFacei];
                rhoVol_[own] -= (rhoSurfPatch[pFacei]*UsurfPatch[pFacei]
                               &SfPatch[pFacei])*dt/V[own];
                rhoUvol[own] -= (rhoSurfPatch[pFacei]*UsurfPatch[pFacei]*UsurfPatch[pFacei]
                                + stressSurfPatch[pFacei])&Sf[pFacei]*dt/V[own];
                scalar rhoEsurf = 
                    rhoSurfPatch[pFacei]
                    *(magSqr(UsurfPatch[pFacei]) + (KInner_ + 3)/2.0*R_.value()*TsurfPatch[pFacei]);

                rhoEvol[own] -= (rhoEsurf*UsurfPatch[pFacei] + qSurfPatch[pFacei])
                               &SfPatch[pFacei]*dt/V[own];
            }
        }
    }

    //- get Prim. from Consv.
    Uvol_ = rhoUvol/rhoVol_;
    Tvol_ = (rhoEvol - 0.5*rhoVol_*magSqr(Uvol_))/((KInner_ + 3)/2.0*R_*rhoVol_);

    //- Correct the macro field boundary conition
    Uvol_.correctBoundaryConditions();
    Tvol_.correctBoundaryConditions();
    //- Note for maxwell wall, the operation here update 
    //- the boundary rho field but it's meaningless.

    //- The vol macro field's boundary field is meanless!
    //rhoVol_.correctBoundaryConditions(); 

    //- update tau
    tauVol_ = updateTau(Tvol_, rhoVol_);

    //- peculiar vel.
    volVectorField c = Uvol_;

    //-get heat flux
    forAll(DV_, dvi)
    {
        discreteVelocity& dv = DV_[dvi];
        c = dv.xi() - Uvol_;
        qVol_ += 0.5*dXiCellSize_*dv.weight()*c
          *(
               magSqr(c)*dv.gTildeVol() 
             + dv.hTildeVol()
           );
    }
    //- correction for bar to original
    qVol_ = 2.0*tauVol_/(2.0*tauVol_ + time_.deltaT()*Pr_)*qVol_;
}


template<template<class> class PatchType, class GeoMesh> 
Foam::tmp<Foam::GeometricField<scalar, PatchType, GeoMesh> >
Foam::fvDVM::updateTau
(
    const GeometricField<scalar, PatchType, GeoMesh>& T, 
    const GeometricField<scalar, PatchType, GeoMesh>& rho
)
{
    tmp<GeometricField<scalar, PatchType, GeoMesh> > tTau
    (
        new GeometricField<scalar, PatchType, GeoMesh>
        (
            IOobject
            (
                "tau",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar( "0", dimTime, 0)
        )
    );
    GeometricField<scalar, PatchType, GeoMesh>& tau = tTau();
    //Debug
    Info << "Tmax = " <<  gMax(T) << endl;
    Info << "Tmin = " <<  gMin(T) << endl;
    tau = muRef_*exp(omega_*log(T/Tref_))/rho/T/R_;
    return tTau;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvDVM::fvDVM
(
    volScalarField& rho,
    volVectorField& U,
    volScalarField& T
)
:
    IOdictionary
    (
        IOobject
        (
            "DVMProperties",
            T.time().constant(),
            T.mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(rho.mesh()),
    time_(rho.time()),
    rhoVol_(rho),
    Uvol_(U),
    Tvol_(T),
    fvDVMparas_(subOrEmptyDict("fvDVMparas")),
    gasProperties_(subOrEmptyDict("gasProperties")),
    nXiPerDim_(readLabel(fvDVMparas_.lookup("nDV"))),
    xiMax_(fvDVMparas_.lookup("xiMax")),
    xiMin_(fvDVMparas_.lookup("xiMin")),
    dXi_((xiMax_-xiMin_)/(nXiPerDim_ - 1)),
    dXiCellSize_
    (
        "dXiCellSize",
        pow(dimLength/dimTime, 3),
        scalar(1.0)
    ),
    macroFlux_(fvDVMparas_.lookupOrDefault("macroFlux", word("no"))),
    res_(fvDVMparas_.lookupOrDefault("res", 1.0e-12)),
    checkSteps_(fvDVMparas_.lookupOrDefault("checkSteps", 100)),
    R_(gasProperties_.lookup("R")),
    omega_(readScalar(gasProperties_.lookup("omega"))),
    Tref_(gasProperties_.lookup("Tref")),
    muRef_(gasProperties_.lookup("muRef")),
    Pr_(readScalar(gasProperties_.lookup("Pr"))),
    KInner_((gasProperties_.lookupOrDefault("KInner", 0))),
    DV_(0),
    rhoSurf_
    (
        IOobject
        (
            "rhoSurf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", rho.dimensions(), 0)
    ),
    Tsurf_
    (
        IOobject
        (
            "Tsurf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", T.dimensions(), 0)
    ),
    Usurf_
    (
        IOobject
        (
            "Usurf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector( "0", U.dimensions(), vector(0,0,0))
    ),
    qSurf_
    (
        IOobject
        (
            "qSurf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector( "0", dimMass/pow(dimTime,3), vector(0,0,0))
    ),
    stressSurf_
    (
        IOobject
        (
            "stressSurf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor("0", dimensionSet(1,-1,-2,0,0,0,0), pTraits<tensor>::zero)
    ),
    qVol_
    (
        IOobject
        (
            "q",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector( "0", dimMass/pow(dimTime,3), vector(0,0,0))
    ),
    tauVol_
    (
        IOobject
        (
            "tauVol",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", dimTime, 0)
    ),
    tauSurf_
    (
        IOobject
        (
            "tauSurf",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar( "0", dimTime, 0)
    )
{
    initialiseDV();
    setCalculatedMaxwellRhoBC();
    updateTau(Tvol_, rhoVol_); //calculate the tau at cell when init
    Usurf_ = fvc::interpolate(Uvol_, "linear"); // for first time Dt calculation.
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvDVM::~fvDVM()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvDVM::evolution()
{
    updateGHbarPvol();
    updateGHbarSurf();
    updateMaxwellWallRho();
    updateGHbarSurfMaxwellWallIn();
    updateMacroSurf();
    updateGHsurf();
    updateGHtildeVol();
    updateMacroVol();
}


void Foam::fvDVM::getCoNum(scalar& maxCoNum, scalar& meanCoNum)
{
    scalar dt = time_.deltaTValue();
    scalarField UbyDx =
        mesh_.surfaceInterpolation::deltaCoeffs()
       *(mag(Usurf_) + sqrt(scalar(mesh_.nSolutionD()))*xiMax_);
    maxCoNum = gMax(UbyDx)*dt;
    meanCoNum = gSum(UbyDx)/UbyDx.size()*dt;
}

// ************************************************************************* //
