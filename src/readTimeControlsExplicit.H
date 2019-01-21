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

Global
    readTimeControls

Author
    lhzhu
Description
    Read the control parameters used by setDeltaT,
    Add a convergence check interval: e.g., 
           convergeCheckSteps 100`
    and a convergence criteria, e.g,
           convergeTol        1.0e-3
    The convergence is check by the temperature field

\*---------------------------------------------------------------------------*/

const bool adjustTimeStep =
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false);

scalar maxCo =
    runTime.controlDict().lookupOrDefault<scalar>("maxCo", 1.0);

scalar maxDeltaT =
    runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

scalar convergeTol =
    runTime.controlDict().lookupOrDefault<scalar>("convergeTol", 1.0e-9);

label convergeCheckSteps =
    runTime.controlDict().lookupOrDefault("convergeCheckSteps", 1);

//TemperatureChange = 1.0;
//rhoChange = 1.0;
//Uchange = 1.0;
// ************************************************************************* //
//
