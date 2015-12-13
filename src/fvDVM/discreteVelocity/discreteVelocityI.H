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

inline const dimensionedVector Foam::discreteVelocity::xi() const
{
    return xi_;
}

inline scalar Foam::discreteVelocity::weight() const
{
    return weight_;
}

inline const volScalarField& Foam::discreteVelocity::gTildeVol() const
{
    return gTildeVol_;
}

inline const volScalarField& Foam::discreteVelocity::hTildeVol() const
{
    return hTildeVol_;
}

inline const volScalarField& Foam::discreteVelocity::gBarPvol() const
{
    return gBarPvol_;
}

inline const volScalarField& Foam::discreteVelocity::hBarPvol() const
{
    return hBarPvol_;
}

inline const surfaceScalarField& Foam::discreteVelocity::gSurf() const
{
    return gSurf_;
}

inline const surfaceScalarField& Foam::discreteVelocity::hSurf() const
{
    return hSurf_;
}


// ************************************************************************* //
