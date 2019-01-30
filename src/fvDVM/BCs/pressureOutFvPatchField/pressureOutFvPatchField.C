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

\*---------------------------------------------------------------------------*/

#include <mpi.h>
#include "pressureOutFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
pressureOutFvPatchField<Type>::pressureOutFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    pressureOut_(0.0)
{}


template<class Type>
pressureOutFvPatchField<Type>::pressureOutFvPatchField
(
    const pressureOutFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    pressureOut_(ptf.pressureOut_)
{}


template<class Type>
pressureOutFvPatchField<Type>::pressureOutFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, false),
    pressureOut_(dict.lookupOrDefault<scalar>("pressureOut", -VGREAT))
{
    //evaluate();
    Field<Type>::operator=(pTraits<Type>::one);
}


template<class Type>
pressureOutFvPatchField<Type>::pressureOutFvPatchField
(
    const pressureOutFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    pressureOut_(ptf.pressureOut_)
{}


template<class Type>
pressureOutFvPatchField<Type>::pressureOutFvPatchField
(
    const pressureOutFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    pressureOut_(ptf.pressureOut_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type>
void pressureOutFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
}


template<class Type>
void pressureOutFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const pressureOutFvPatchField<Type>& fgptf =
        refCast<const pressureOutFvPatchField<Type> >(ptf);

    pressureOut_ = fgptf.pressureOut_;
}

template<class Type>
void pressureOutFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
    
    // Do nothing

    //Info << inComingByRho_ << endl;
    fvPatchField<Type>::evaluate();
}


template<class Type>
void pressureOutFvPatchField<Type>::write(Ostream& os) const
{
        fvPatchField<Type>::write(os);
        this->writeEntry("value", os);
        os.writeKeyword("pressureOut")
            << pressureOut_ << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
