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
    Foam::symmetryModFvPatchField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition supplies a fixed gradient condition, such that
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
        type            symmetryMod;
        gradient        uniform 0;
    }
    \endverbatim

SourceFiles
    symmetryModFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef symmetryModFvPatchField_H
#define symmetryModFvPatchField_H

#include "foam_defs.h"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class symmetryModFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>

class symmetryModFvPatchField
:
    public fvPatchField<Type>
{
    // Private data

        Field<Type> dfContainer_;

public:

    //- Runtime type information
    TypeName("symmetryMod");


    // Constructors

        //- Construct from patch and internal field
        symmetryModFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        symmetryModFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given symmetryModFvPatchField
        //  onto a new patch
        symmetryModFvPatchField
        (
            const symmetryModFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        symmetryModFvPatchField
        (
            const symmetryModFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new symmetryModFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        symmetryModFvPatchField
        (
            const symmetryModFvPatchField<Type>&,
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
                new symmetryModFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        // Return defining fields

            //- Return inComingByRho at boundary
            virtual Field<Type>& dfContainer()
            {
                return dfContainer_;
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
#   include "symmetryModFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
