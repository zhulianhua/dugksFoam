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

Class
    Foam::pressureInFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition supplies a fixed pressure condition, such that
    the patch values are calculated using:

        \f[
            x_p = x_c + \frac{\nabla(x)}{\Delta}
        \f]

    where
    \vartable
        x_p      | patch values
        x_c      | internal field values
        \nabla(x)| gradient (user-specified)
        \Delta   | inverse distance from patch face centre to cell centre
    \endvartable

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        gradient     | gradient                | yes         | 
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            pressureIn;
        gradient        uniform 0;
    }
    \endverbatim

SourceFiles
    pressureInFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef pressureInFvPatchField_H
#define pressureInFvPatchField_H

#include "foam_defs.h"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class pressureInFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>

class pressureInFvPatchField
:
    public fvPatchField<Type>
{
    // Private data

        scalar pressureIn_;

public:

    //- Runtime type information
    TypeName("pressureIn");


    // Constructors

        //- Construct from patch and internal field
        pressureInFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pressureInFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given pressureInFvPatchField
        //  onto a new patch
        pressureInFvPatchField
        (
            const pressureInFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        pressureInFvPatchField
        (
            const pressureInFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new pressureInFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        pressureInFvPatchField
        (
            const pressureInFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new pressureInFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Return defining fields

            //- Return inComingByRho at boundary
            virtual scalar pressureIn()
            {
                return pressureIn_;
            }

            //- Evaluate the patch field
            virtual void evaluate
            (
#if FOAM_MAJOR <= 4
                const Pstream::commsTypes commsType=Pstream::blocking
#else
                const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
#endif
            );

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchField<Type>&,
                const labelList&
            );

        //- Write
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "pressureInFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
