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
#include "fvDVM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readTimeControlsExplicit.H"

    fvDVM dvm(rho, U, T);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    label It = 0;
    while (runTime.run() && (TemperatureChange > convergeTol))
    {
        #include "CourantNo.H" // calculate the Co num
        #include "readTimeControlsExplicit.H"
        #include "setDeltaT.H"

        runTime++;
        It++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        dvm.evolution();

        runTime.write();

        Info<< "Step =" << It << "  ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        if(It%convergeCheckSteps == 0 && It >= convergeCheckSteps)
        {
            tmp<Foam::GeometricField<scalar, Foam::fvPatchField, Foam::volMesh> > 
                deltaTem = mag(T-Told);
            TemperatureChange = gSum(deltaTem())/gSum(T);
            Info << "Temperature changes = " << TemperatureChange << nl << endl;
            Told = T;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
