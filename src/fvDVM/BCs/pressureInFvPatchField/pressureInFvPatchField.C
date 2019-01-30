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
#include "pressureInFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
pressureInFvPatchField<Type>::pressureInFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    pressureIn_(0.0)
{}


template<class Type>
pressureInFvPatchField<Type>::pressureInFvPatchField
(
    const pressureInFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    pressureIn_(ptf.pressureIn_)
{}


template<class Type>
pressureInFvPatchField<Type>::pressureInFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, false),
    pressureIn_(dict.lookupOrDefault<scalar>("pressureIn", -VGREAT))
{
    //evaluate();
    Field<Type>::operator=(pTraits<Type>::one);
}


template<class Type>
pressureInFvPatchField<Type>::pressureInFvPatchField
(
    const pressureInFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    pressureIn_(ptf.pressureIn_)
{}


template<class Type>
pressureInFvPatchField<Type>::pressureInFvPatchField
(
    const pressureInFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    pressureIn_(ptf.pressureIn_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void pressureInFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
}


template<class Type>
void pressureInFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const pressureInFvPatchField<Type>& fgptf =
        refCast<const pressureInFvPatchField<Type> >(ptf);

    pressureIn_ = fgptf.pressureIn_;
}

template<class Type>
void pressureInFvPatchField<Type>::evaluate(const Pstream::commsTypes)
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
void pressureInFvPatchField<Type>::write(Ostream& os) const
{
        fvPatchField<Type>::write(os);
        this->writeEntry("value", os);
        os.writeKeyword("pressureIn")
            << pressureIn_ << token::END_STATEMENT << nl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
