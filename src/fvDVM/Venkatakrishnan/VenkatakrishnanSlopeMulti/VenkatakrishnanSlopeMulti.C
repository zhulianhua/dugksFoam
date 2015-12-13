/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::VenkatakrishnanSlopeMultiLimiter

Description
    Class with limiter function which returns the limiter
    VenkatakrishnanSlopeMulti
    V. Venkatakrishnan. "Convergence to Steady State Solutions of the Euler
    Equations on Unstructured Grids with Limiters."
    Journal of Computational Physics, 118:120–130, 1995.

    Blazek, Jiri. "Computational Fluid Dynamics: Principles and Applications",
    2001, Referex Engineering, Elsevier

Author
    Oliver Borm  All rights reserved.

SourceFiles
    VenkatakrishnanSlopeMulti.H

\*---------------------------------------------------------------------------*/

#ifndef VenkatakrishnanSlopeMulti_H
#define VenkatakrishnanSlopeMulti_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class VenkatakrishnanSlopeMultiLimiter Declaration
\*---------------------------------------------------------------------------*/

class VenkatakrishnanSlopeMultiLimiter
{
    scalar k_;

public:

    VenkatakrishnanSlopeMultiLimiter(Istream& is)
    :
        k_(readScalar(is))
    {
        if (k_ < 0 )
        {
            FatalIOErrorIn("VenkatakrishnanSlopeMultiLimiter(Istream& is)", is)
                << "coefficient = " << k_
                << " should be >= 0"
                << exit(FatalIOError);
        }
    }

    inline scalar k() const 
    {
        return k_;
    }
    // Venkatakrishnan Slope Limiter
    // Implementation acccording to Blazek Chap.5 / p. 167
    inline scalar limitFace
    (
        const scalar cellVolume,
        const scalar deltaOneMax,
        const scalar deltaOneMin,
        const scalar deltaTwo
    ) const
    {
        // is always > 0 if k_ is > 0
        scalar sqrEpsilon = pow3(k_)*cellVolume;
        scalar twoSqrDeltaTwo = 2.0*sqr(deltaTwo);
        scalar sqrDeltaOneMax = sqr(deltaOneMax);
        scalar sqrDeltaOneMin = sqr(deltaOneMin);

        // deltaOneMax = phiMax - phiP is always >= 0
        // deltaOneMax = phiMin - phiP is always <= 0
        // deltaTwo can be positive or negative,
        // if deltaTwo > 0, then scalar positive should be always > 0
        // if deltaTwo < 0, then scalar negative should be always < 0

//         scalar positive = max(((sqrDeltaOneMax+sqrEpsilon)*deltaTwo+twoSqrDeltaTwo*deltaOneMax)
//             /stabilise((sqrDeltaOneMax+twoSqrDeltaTwo+deltaOneMax*deltaTwo+sqrEpsilon),VSMALL),0);
// 
//         scalar negative = min(((sqrDeltaOneMin+sqrEpsilon)*deltaTwo+twoSqrDeltaTwo*deltaOneMin)
//             /stabilise((sqrDeltaOneMin+twoSqrDeltaTwo+deltaOneMin*deltaTwo+sqrEpsilon),VSMALL),0);
// 
//         deltaTwo = stabilise(deltaTwo,VSMALL);
// 
//         // This sets the limiter to zero if the gradient is zero
//         // and not to One as proposed by the author. But if the gradient is
//         // already zero the slope limiter could be arbitrary
//         return min(max(positive/deltaTwo,negative/deltaTwo),1);

////////////////////////////

        if (deltaTwo > 0.0)
        {
            return ((sqrDeltaOneMax+sqrEpsilon)+2.0*deltaTwo*deltaOneMax)
            /stabilise((sqrDeltaOneMax+twoSqrDeltaTwo+deltaOneMax*deltaTwo+sqrEpsilon),VSMALL);
        }
        else if (deltaTwo < 0.0)
        {
            return ((sqrDeltaOneMin+sqrEpsilon)+2.0*deltaTwo*deltaOneMin)
            /stabilise((sqrDeltaOneMin+twoSqrDeltaTwo+deltaOneMin*deltaTwo+sqrEpsilon),VSMALL);
        }
        else
        {
            return 1.0;
        }
//         deltaTwo = stabilise(deltaTwo,VSMALL);

        // This sets the limiter to zero if the gradient is zero
        // and not to One as proposed by the author. But if the gradient is
        // already zero the slope limiter could be arbitrary
//         return min(max(positive,negative),1);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
