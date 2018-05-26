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

#include <mpi.h>
#include "fvDVM.H"
#include "constants.H"
#include "fvm.H"
#include "calculatedMaxwellFvPatchField.H"
#include "symmetryModFvPatchField.H"
#include "pressureInFvPatchField.H"
#include "pressureOutFvPatchField.H"
#include "scalarIOList.H"
#include "fieldMPIreducer.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

#if FOAM_MAJOR <= 3
    #define BOUNDARY_FIELD_REF boundaryField()
#else
    #define BOUNDARY_FIELD_REF boundaryFieldRef()
#endif


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

    scalarField weightsGlobal;
    vectorField XisGlobal;
    labelField  symmXtgID;
    labelField  symmYtgID;
    labelField  symmZtgID;

    if (mesh_.nSolutionD() == 3)    //3D(X & Y & Z)
    {
        nXiX_ = nXiY_ = nXiZ_ = nXiPerDim_;
        nXi_ = nXiX_*nXiY_*nXiZ_;

        weightsGlobal.setSize(nXi_);
        XisGlobal.setSize(nXi_);
        symmXtgID.setSize(nXi_);
        symmYtgID.setSize(nXi_);
        symmZtgID.setSize(nXi_);

        label i = 0;
        for (label iz = 0; iz < nXiZ_; iz++)
        {
            for (label iy = 0; iy < nXiY_; iy++)
            {
                for (label ix = 0; ix < nXiZ_; ix++)
                {
                    scalar weight = weights1D[iz]*weights1D[iy]*weights1D[ix];
                    vector xi(Xis[ix], Xis[iy], Xis[iz]);
                    weightsGlobal[i] = weight;
                    XisGlobal[i] = xi;
                    symmXtgID[i] = iz*nXiY_*nXiX_ + iy*nXiX_ + (nXiX_ - ix -1);
                    symmYtgID[i] = iz*nXiY_*nXiX_ + (nXiY_ - iy - 1)*nXiX_ + ix;
                    symmZtgID[i] = (nXiZ_ - iz -1)*nXiY_*nXiX_ + iy*nXiX_ + ix;
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
            weightsGlobal.setSize(nXi_);
            XisGlobal.setSize(nXi_);
            symmXtgID.setSize(nXi_);
            symmYtgID.setSize(nXi_);
            symmZtgID.setSize(nXi_);
            label i = 0;
            for (label iy = 0; iy < nXiY_; iy++)
            {
                for (label ix = 0; ix < nXiX_; ix++)
                {
                    scalar weight = weights1D[iy]*weights1D[ix]*1;
                    vector xi(Xis[ix], Xis[iy], 0.0);
                    weightsGlobal[i] = weight;
                    XisGlobal[i] = xi;
                    symmXtgID[i] = iy*nXiX_ + (nXiX_ - ix -1);
                    symmYtgID[i] = (nXiY_ - iy - 1)*nXiX_ + ix;
                    symmZtgID[i] = 0;
                    i++;
                }
            }
        }
        else    //1D (X)
        {
            nXiX_ = nXiPerDim_;
            nXiY_ = nXiZ_ = 1;
            nXi_ = nXiX_*nXiY_*nXiZ_;
            weightsGlobal.setSize(nXi_);
            XisGlobal.setSize(nXi_);
            symmXtgID.setSize(nXi_);
            symmYtgID.setSize(nXi_);
            symmZtgID.setSize(nXi_);
            label i = 0;
            for (label ix = 0; ix < nXiX_; ix++)
            {
                scalar weight = weights1D[ix]*1*1;
                vector xi(Xis[ix], 0.0, 0.0);
                weightsGlobal[i] = weight;
                XisGlobal[i] = xi;
                symmXtgID[i] = (nXiX_ - ix -1);
                symmYtgID[i] = 0;              
                symmZtgID[i] = 0;              
                i++;
            }
        }
    }

    if(mpiReducer_.rank() == 0)
    {
        Info<< "fvDVM : Allocated " << XisGlobal.size()
            << " discrete velocities" << endl;
    }
    label nA = nXi_ / mpiReducer_.nproc();
    label nB = nXi_ - nA*mpiReducer_.nproc();
    label nXiPart = nA + (label)(mpiReducer_.rank() < nB);
    DV_.setSize(nXiPart);
    if(mpiReducer_.rank() == 0)
    {
        Info << "nproc    " << mpiReducer_.nproc()  << endl;
        Info << "nXisPart " << nXiPart << endl;
    }
    label chunk = 0;
    label gid = 0;

    forAll(DV_, i)
    {
        gid = chunk + mpiReducer_.rank();
        DV_.set
        (
            i,
            new discreteVelocity
            (
                *this,
                mesh_,
                time_,
                weightsGlobal[gid],
                dimensionedVector( "xi", dimLength/dimTime, XisGlobal[gid]),
                i,
                symmXtgID[gid],
                symmYtgID[gid],
                symmZtgID[gid]
            )
        );
        chunk += mpiReducer_.nproc();
    }
}

void Foam::fvDVM::setCalculatedMaxwellRhoBC()
{
#if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
#else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
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

            if(args_.optionFound("dvParallel"))
                mpiReducer_.reduceField(rhoPatch.inComingByRho());
        }

    }
}

void Foam::fvDVM::setSymmetryModRhoBC()
{
    //prepare the container (set size) to store all DF on the patchi
#if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
#else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
    forAll(rhoBCs, patchi)
    {
        label ps = rhoBCs[patchi].size();
        if (rhoBCs[patchi].type() == "symmetryMod")
        {
            symmetryModFvPatchField<scalar>& rhoPatch = 
                refCast<symmetryModFvPatchField<scalar> >(rhoBCs[patchi]);
            rhoPatch.dfContainer().setSize(ps*nXi_*2); //*2 means g and h
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
#if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
#else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
    forAll(rhoBCs, patchi)
    {
        if (rhoBCs[patchi].type() == "calculatedMaxwell")
        {
            calculatedMaxwellFvPatchField<scalar>& rhoPatch = 
                refCast<calculatedMaxwellFvPatchField<scalar> >(rhoBCs[patchi]);
            if(args_.optionFound("dvParallel"))
                mpiReducer_.reduceField(rhoPatch.outGoing());
        }
    }
    rhoVol_.correctBoundaryConditions();
}

void Foam::fvDVM::updateGHbarSurfMaxwellWallIn()
{
    forAll(DV_, DVid)
        DV_[DVid].updateGHbarSurfMaxwellWallIn();
}

void Foam::fvDVM::updateGHbarSurfSymmetryIn()
{
    //1. copy all DV's g/h to rho patch's dfContainer
    //2. MPI_Allgather the rho patch's dfContainer
    //if(args_.optionFound("dvParallel"))
    //{
    label rank  = mpiReducer_.rank();
    label nproc = mpiReducer_.nproc();
#if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
#else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
    forAll(rhoBCs, patchi)
    {
        label ps = rhoBCs[patchi].size();
        if (rhoBCs[patchi].type() == "symmetryMod")
        {
            symmetryModFvPatchField<scalar>& rhoPatch = 
                refCast<symmetryModFvPatchField<scalar> >(rhoBCs[patchi]);
            //compose the recvcout and displacement array
            labelField recvc(nproc);
            labelField displ(nproc);
            label chunck = nXi_/nproc;
            label left   = nXi_%nproc;
            forAll(recvc, i)
            {
                recvc[i] = 2*ps*(chunck + (i<left)) ;
                if(i<=left)
                    displ[i] = i*2*ps*(chunck + 1); // (i<=nXi_%nproc)
                else
                    displ[i] = 2*ps*(left*(chunck +1) + (i-left)*(chunck));
            }

            // check 12*28+15 dv's g
            label did = 1709;
            label pp  = did%nproc;
            label lid = did/nproc;
            //if(rank==pp)
            //{
                //Info << "processing by rank " << rank << endl;
                //Info << "12*28+15 outging g " << DV_[lid].gSurf()[0] << endl;
                //Info << "12*28+15 outging xi " << DV_[lid].xi() <<endl;
                //Info << "12*28+15 outging at boundary " << DV_[lid].gSurf().boundaryField()[patchi][0] << endl;
            //}
            // memcpy each dv's g/h to rho
            forAll(DV_, DVid)
            {
                //label shift = (nXi_ / nproc * rank + DVid)*2*ps;
                label shift = displ[rank] + DVid*2*ps;
                memcpy( (rhoPatch.dfContainer().data() + shift),
                        DV_[DVid].gSurf().boundaryField()[patchi].cdata(), ps*sizeof(scalar));
                memcpy( (rhoPatch.dfContainer().data() + shift + ps),
                        DV_[DVid].hSurf().boundaryField()[patchi].cdata(), ps*sizeof(scalar));
            }

            // check 
            //if(rank == pp)
                //Info << "dv gid 1709's g = " <<rhoPatch.dfContainer()[displ[pp]+lid*2*ps+32]<< endl;;


            //Allgather
            MPI_Allgatherv(
                //rhoPatch.dfContainer().data() + displ[rank],//2*ps*nXI_/nproc*rank, //send*
                MPI_IN_PLACE,
                2*ps*DV_.size(), //(how many DV i processed) * 2 * patch size
                MPI_DOUBLE,
                rhoPatch.dfContainer().data(),
                recvc.data(),
                displ.data(),
                MPI_DOUBLE,
                MPI_COMM_WORLD
                );
        }
    }
        forAll(DV_, DVid)
            DV_[DVid].updateGHbarSurfSymmetryIn();
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
    //if(Pstream::parRun())
    //{
        //Pout << "My proc NO " << Pstream::myProcNo() << endl;
        //if(Pstream::myProcNo() == 1)
            //Pout << "g[78] at proc 1  interface" <<  DV_[100].gSurf().boundaryField()[mesh_.boundaryMesh().findPatchID("procBoundary1to0")][0] << endl;
        //if(Pstream::myProcNo() == 0)
            //Pout << "g[78] at proc 0  interface" <<  DV_[100].gSurf().boundaryField()[mesh_.boundaryMesh().findPatchID("procBoundary0to1")][0] << endl;
    //}else
    //{
        //Info << "g[78] at 50 th  interface" <<  DV_[100].gSurf()[49] << endl;
    //}

    // Conserved variable, now zero as the prime variable 
    // has been set to zero
    surfaceVectorField rhoUsurf = rhoSurf_*Usurf_;
    surfaceScalarField rhoEsurf = rhoSurf_*magSqr(Usurf_);

    forAll(DV_, dvi)
    {
        discreteVelocity& dv = DV_[dvi];
        rhoSurf_  += dXiCellSize_*dv.weight()*dv.gSurf();
        rhoUsurf  += dXiCellSize_*dv.weight()*dv.gSurf()*dv.xi();
        rhoEsurf  += 0.5*dXiCellSize_*dv.weight()
           *(
                dv.gSurf()*magSqr(dv.xi()) 
              + dv.hSurf()
            );
    }

    if(args_.optionFound("dvParallel"))
    {
        mpiReducer_.reduceField(rhoSurf_);
        mpiReducer_.reduceField(rhoUsurf);
        mpiReducer_.reduceField(rhoEsurf);
    }

    //- get Prim. from Consv.
    Usurf_ = rhoUsurf/rhoSurf_;

    Tsurf_ = (rhoEsurf - 0.5*rhoSurf_*magSqr(Usurf_))/((KInner_ + 3)/2.0*R_*rhoSurf_);

    tauSurf_ = updateTau(Tsurf_, rhoSurf_);
    //- peculiar vel.

    surfaceVectorField c = Usurf_;

    //-get part heat flux 
    forAll(DV_, dvi)
    {
        discreteVelocity& dv = DV_[dvi];
        c = dv.xi() - Usurf_;
        qSurf_ += 0.5*dXiCellSize_*dv.weight()*c
            *(
                 magSqr(c)*dv.gSurf() 
               + dv.hSurf()
             );
        //- stressSurf is useless as we never update cell macro by macro flux 
        //- Comment out it as it is expansive
        //stressSurf_ += 
            //dXiCellSize_*dv.weight()*dv.gSurf()*c*c;
    }
    //- Get global heat flux, via MPI_Allreuce
    if(args_.optionFound("dvParallel"))
        mpiReducer_.reduceField(qSurf_);

    //- correction for bar to original
    qSurf_ = 2.0*tauSurf_/(2.0*tauSurf_ + 0.5*time_.deltaT()*Pr_)*qSurf_;



    //- stress at surf is not used, as we dont't update macro in cell by macro flux at surface
    //stressSurf_ = 
        //2.0*tauSurf_/(2.0*tauSurf_ + 0.5*time_.deltaT())*stressSurf_;

    //- heat flux at wall is specially defined. as it ignores the velocity and temperature slip
    //- NOTE: To be changed as it is part macro, but it will not affect the innner fields, so we change it later
#if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
#else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
    qWall_ = dimensionedVector("0", qWall_.dimensions(), vector(0, 0, 0));
    stressWall_ = dimensionedTensor
        (
            "0", 
            stressWall_.dimensions(), 
            pTraits<tensor>::zero
        );
    forAll(rhoBCs, patchi)
    {
        if (rhoBCs[patchi].type() == "calculatedMaxwell")
        {
            fvPatchField<vector>& qPatch = qWall_.BOUNDARY_FIELD_REF[patchi];
            fvPatchField<vector>& Upatch = Uvol_.BOUNDARY_FIELD_REF[patchi];
            fvPatchField<tensor>& stressPatch = stressWall_.BOUNDARY_FIELD_REF[patchi];
            //- tau at surface use the tau at slip temperature as it is.
            fvsPatchField<scalar>&  tauPatch = tauSurf_.BOUNDARY_FIELD_REF[patchi];
            forAll(qPatch, facei)
            {
                forAll(DV_, dvi)
                {
                    scalar dXiCellSize = dXiCellSize_.value();
                    discreteVelocity& dv = DV_[dvi];
                    vector xi = dv.xi().value();
                    vector c = xi - Upatch[facei];
                    qPatch[facei] += 0.5*dXiCellSize*dv.weight()*c  //sometimes wall moves, then c != \xi
                        *(
                             magSqr(c)*dv.gSurf().boundaryField()[patchi][facei]
                           + dv.hSurf().boundaryField()[patchi][facei]
                         );
                    stressPatch[facei] += 
                        dXiCellSize*dv.weight()*dv.gSurf().boundaryField()[patchi][facei]*xi*xi;
                }
                qPatch[facei] = 2.0*tauPatch[facei]/(2.0*tauPatch[facei] + 0.5*time_.deltaT().value()*Pr_)*qPatch[facei];
                stressPatch[facei] = 
                    2.0*tauPatch[facei]/(2.0*tauPatch[facei] + 0.5*time_.deltaT().value())*stressPatch[facei];
            }
            if(args_.optionFound("dvParallel"))
            {
                mpiReducer_.reduceField(qPatch);
                mpiReducer_.reduceField(stressPatch);
            }
        }
    }
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
    //- Old macros, used only if we update using macro fluxes.
    volVectorField rhoUvol = rhoVol_*Uvol_;
    volScalarField rhoEvol = rhoVol_*(magSqr(Uvol_) + (KInner_ + 3)/2.0*R_*Tvol_);
    qVol_ = dimensionedVector("0", qVol_.dimensions(), vector(0, 0, 0));

    if(macroFlux_ == "no") // update cell macro by moment from DF
    {
        //- init to zeros
        rhoVol_ = dimensionedScalar("0", rhoVol_.dimensions(), 0);
        rhoUvol = dimensionedVector("0", rhoUvol.dimensions(), vector(0, 0, 0));
        rhoEvol = dimensionedScalar("0", rhoEvol.dimensions(), 0);

        //- get part macro
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
        //- get global macro via MPI_Allreduce
        if(args_.optionFound("dvParallel"))
        {
            mpiReducer_.reduceField(rhoVol_);
            mpiReducer_.reduceField(rhoUvol);
            mpiReducer_.reduceField(rhoEvol);
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

    //-get part heat flux
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

    //- get global heat flux via MPI_Allreduce
    if(args_.optionFound("dvParallel"))
        mpiReducer_.reduceField(qVol_);
    //- correction for bar to original
    qVol_ = 2.0*tauVol_/(2.0*tauVol_ + time_.deltaT()*Pr_)*qVol_;
}

void Foam::fvDVM::updatePressureInOutBC()
{
    // for pressureIn and pressureOut BC, the boundary value of Uvol(in/out) and Tvol(in/out) should be updated here!
    // boundary faces
#if FOAM_MAJOR <= 3
    GeometricField<scalar, fvPatchField, volMesh>::GeometricBoundaryField& 
        rhoBCs = rhoVol_.boundaryField();
#else
    GeometricField<scalar, fvPatchField, volMesh>::Boundary& 
        rhoBCs = rhoVol_.BOUNDARY_FIELD_REF;
#endif
    forAll(rhoBCs, patchi)
    {
        if (rhoBCs[patchi].type() == "pressureIn")
        {
            const fvsPatchField<vector>& SfPatch = mesh_.Sf().boundaryField()[patchi];
            const fvsPatchField<scalar>& magSfPatch = mesh_.magSf().boundaryField()[patchi];
            pressureInFvPatchField<scalar>& rhoPatch = 
                refCast<pressureInFvPatchField<scalar> >(rhoBCs[patchi]);
            fvPatchField<vector>& Upatch = Uvol_.BOUNDARY_FIELD_REF[patchi];
            const fvPatchField<scalar>& Tpatch = Tvol_.boundaryField()[patchi];
            const scalar pressureIn = rhoPatch.pressureIn();
            // now changed rho and U patch
            const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
            forAll(rhoPatch, facei)
            {
                const scalar  Tin = Tpatch[facei];
                // change density
                rhoPatch[facei] = pressureIn/R_.value()/Tin; // Accturally not changed at all :p

                // inner boundary cell data state data
                label own = pOwner[facei];
                vector Ui = Uvol_[own];
                scalar Ti = Tvol_[own];
                scalar rhoi = rhoVol_[own];
                scalar ai = sqrt(R_.value() * Ti * (KInner_ + 5)/(KInner_ + 3)); // sos

                // change normal velocity component based on the characteristics
                vector norm = SfPatch[facei]/magSfPatch[facei]; // boundary face normal vector
                scalar Un = Ui & norm; // normal component
                scalar UnIn = Un + (pressureIn - rhoi * R_.value() * Ti)/rhoi/ai; // change normal component
                Upatch[facei] = UnIn * norm + (Ui - Un * norm); // tangential component not changed.
            }
        }
        else if(rhoBCs[patchi].type() == "pressureOut")
        {
            const fvsPatchField<vector>& SfPatch = mesh_.Sf().boundaryField()[patchi];
            const fvsPatchField<scalar>& magSfPatch = mesh_.magSf().boundaryField()[patchi];
            pressureOutFvPatchField<scalar>& rhoPatch = 
                refCast<pressureOutFvPatchField<scalar> >(rhoBCs[patchi]);
            fvPatchField<vector>& Upatch = Uvol_.BOUNDARY_FIELD_REF[patchi];
            fvPatchField<scalar>& Tpatch = Tvol_.BOUNDARY_FIELD_REF[patchi];
            const scalar pressureOut = rhoPatch.pressureOut();
            // now changed rho and U patch
            const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();
            forAll(rhoPatch, facei)
            {
                // inner cell data state data
                label own = pOwner[facei];
                vector Ui = Uvol_[own];
                scalar Ti = Tvol_[own];
                scalar rhoi = rhoVol_[own];
                scalar ai = sqrt(R_.value() * Ti * (KInner_ + 5)/(KInner_ + 3)); // sos

                // change outlet density
                rhoPatch[facei] = rhoi  +  (pressureOut - rhoi * R_.value() * Ti)/ai/ai; // Accturally not changed at all :p
                Tpatch[facei] = pressureOut/(R_.value() * rhoi);

                // change normal velocity component based on the characteristics
                vector norm = SfPatch[facei]/magSfPatch[facei]; // boundary face normal vector
                scalar Un = Ui & norm; // normal component
                scalar UnIn = Un + ( rhoi * R_.value() * Ti - pressureOut)/rhoi/ai; // change normal component
                Upatch[facei] = UnIn * norm + (Ui - Un * norm); // tangential component not changed.
            }
        }
    }
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
    GeometricField<scalar, PatchType, GeoMesh> tau = tTau();
    tau = muRef_*exp(omega_*log(T/Tref_))/rho/T/R_;
    return tTau;
}


void Foam::fvDVM::writeDF(label cellId)
{
    std::ostringstream convert;
    convert << cellId;
    scalarIOList df
    (
        IOobject
        (
             "DF"+convert.str(),
             time_.time().timeName(),
             mesh_,
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
        )
    );
    //set size of df
    df.setSize(nXi_);

    scalarList dfPart(DV_.size());
    //put cellId's DF to dfPart
    forAll(dfPart, dfi)
        dfPart[dfi] = DV_[dfi].gTildeVol()[cellId];

    label nproc = mpiReducer_.nproc();
    //gather
    //tmp list for recv
    scalarList dfRcv(nXi_);

    //Compose displc and recvc
    labelField recvc(nproc);
    labelField displ(nproc);
    label chunck = nXi_/nproc;
    label left   = nXi_%nproc;
    forAll(recvc, i)
    {
        recvc[i] = chunck + (i<left) ;
        if(i<=left)
            displ[i] = i*(chunck + 1); // (i<=nXi_%nproc)
        else
            displ[i] = left*(chunck +1) + (i-left)*(chunck);
    }
    MPI_Gatherv(dfPart.cdata(), dfPart.size(), MPI_DOUBLE,
            dfRcv.data(), recvc.cdata(), displ.cdata(), 
            MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //reposition
    if(mpiReducer_.rank() == 0)
    {
        forAll(df, i)
        {
            label p   = i%nproc;
            label ldi = i/nproc;
            df[i] = dfRcv[displ[p]+ldi];
        }
        df.write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    Foam::fvDVM::fvDVM
(
 volScalarField& rho,
 volVectorField& U,
 volScalarField& T,
 int* argc,
 char*** argv,
 Foam::argList& args
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
        args_(args),
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
        //res_(fvDVMparas_.lookupOrDefault("res", 1.0e-12)),
        //checkSteps_(fvDVMparas_.lookupOrDefault("checkSteps", 100)),
        R_(gasProperties_.lookup("R")),
        omega_(readScalar(gasProperties_.lookup("omega"))),
        Tref_(gasProperties_.lookup("Tref")),
        muRef_(gasProperties_.lookup("muRef")),
        Pr_(readScalar(gasProperties_.lookup("Pr"))),
        KInner_((gasProperties_.lookupOrDefault("KInner", 0))),
        mpiReducer_(args, argc, argv), // args comes from setRootCase.H in dugksFoam.C;
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
          IOobject::AUTO_WRITE
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
          IOobject::AUTO_WRITE
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
        ),
        qWall_
        (
         IOobject
         (
          "qWall",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         mesh_,
         dimensionedVector( "0", dimMass/pow(dimTime,3), vector(0,0,0))
        ),
        stressWall_
        (
         IOobject
         (
          "stressWall",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
         ),
         mesh_,
         dimensionedTensor("0", dimensionSet(1,-1,-2,0,0,0,0), pTraits<tensor>::zero)
        )
{
    initialiseDV();
    setCalculatedMaxwellRhoBC();
    setSymmetryModRhoBC();
    // set initial rho in pressureIn/Out BC
    updatePressureInOutBC();
    tauVol_ = updateTau(Tvol_, rhoVol_); //calculate the tau at cell when init
    Usurf_ = fvc::interpolate(Uvol_, "linear"); // for first time Dt calculation.
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvDVM::~fvDVM()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvDVM::evolution()
{
    //Info << "Begin evolution" << endl;
    updateGHbarPvol();
    //Info << "Done updateGHbarPvol " << endl;
    updateGHbarSurf();
    //Info << "Done updateGHbarSurf " << endl;
    updateMaxwellWallRho();
    //Info << "Done updateMaxwellWallRho " << endl;
    updateGHbarSurfMaxwellWallIn();
    //Info << "Done updateGHbarSurfMaxwellWallIn " << endl;
    updateGHbarSurfSymmetryIn();
    //Info << "Done updateGHbarSurfSymmetryIn " << endl;
    updateMacroSurf();
    //Info << "Done updateMacroSurf " << endl;
    updateGHsurf();
    //Info << "Done updateGHsurf " << endl;
    updateGHtildeVol();
    //Info << "Done updateGHtildeVol " << endl;
    updateMacroVol();
    //
    updatePressureInOutBC();
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


const fieldMPIreducer& Foam::fvDVM::mpiReducer() const
{
    return mpiReducer_;
}


// ************************************************************************* //
