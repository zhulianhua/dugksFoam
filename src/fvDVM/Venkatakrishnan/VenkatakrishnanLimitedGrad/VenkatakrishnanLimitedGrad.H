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

Class
    Foam::fv::VenkatakrishnanLimitedGrad

Description
    VenkatakrishnanLimitedGrad gradient scheme applied to a runTime selected base gradient
    scheme.

    The scalar limiter based on limiting the extrapolated face values
    between the maximum and minumum cell and cell neighbour values and is
    applied to all components of the gradient.

SourceFiles
    VenkatakrishnanLimitedGrad.C

\*---------------------------------------------------------------------------*/

#ifndef VenkatakrishnanLimitedGrad_H
#define VenkatakrishnanLimitedGrad_H

#include "gradScheme.H"
#include "VenkatakrishnanSlopeMulti.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

/*---------------------------------------------------------------------------*\
                       Class VenkatakrishnanLimitedGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class VenkatakrishnanLimitedGrad
:
    public fv::gradScheme<Type>
{
    // Private Data



        //- Limiter coefficient
        tmp<fv::gradScheme<Type> > basicGradScheme_;

        //- Venkatakrishnan slope for every edge
        VenkatakrishnanSlopeMultiLimiter VenSlope_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        VenkatakrishnanLimitedGrad(const VenkatakrishnanLimitedGrad&);

        //- Disallow default bitwise assignment
        void operator=(const VenkatakrishnanLimitedGrad&);


public:

    //- RunTime type information
    TypeName("VenkatakrishnanLimited");


    // Constructors

        //- Construct from mesh and schemeData
        VenkatakrishnanLimitedGrad(const fvMesh& mesh, Istream& schemeData)
        :
            gradScheme<Type>(mesh),
            basicGradScheme_(fv::gradScheme<Type>::New(mesh, schemeData)),
            VenSlope_(schemeData)
        {
        }


    // Member Functions

        //- Return the gradient of the given field to the gradScheme::grad
        //  for optional caching
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
        > calcGrad
        (
            const GeometricField<Type, fvPatchField, volMesh>& vsf,
            const word& name
        ) const;
};


// * * * * * * * * * * * * Inline Member Function  * * * * * * * * * * * * * //



// * * * * * * * * Template Member Function Specialisations  * * * * * * * * //

template<>
tmp<volVectorField> VenkatakrishnanLimitedGrad<scalar>::calcGrad
(
    const volScalarField& vsf,
    const word& name
) const;


template<>
tmp<volTensorField> VenkatakrishnanLimitedGrad<vector>::calcGrad
(
    const volVectorField& vsf,
    const word& name
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
