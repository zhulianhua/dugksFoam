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
#include "calculatedMaxwellFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
calculatedMaxwellFvPatchField<Type>::calculatedMaxwellFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    inComingByRho_(p.size(), pTraits<Type>::zero),
    outGoing_(p.size(), pTraits<Type>::zero)
{}


template<class Type>
calculatedMaxwellFvPatchField<Type>::calculatedMaxwellFvPatchField
(
    const calculatedMaxwellFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    inComingByRho_(ptf.inComingByRho_, mapper),
    outGoing_(ptf.inComingByRho_, mapper)
{}


template<class Type>
calculatedMaxwellFvPatchField<Type>::calculatedMaxwellFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    inComingByRho_(p.size(), pTraits<Type>::zero),
    outGoing_(p.size(), pTraits<Type>::zero)
{
    //evaluate();
    //set field to scalar(1)
    //*this == scalar(1.0);
    Field<Type>::operator=(pTraits<Type>::one);
}


template<class Type>
calculatedMaxwellFvPatchField<Type>::calculatedMaxwellFvPatchField
(
    const calculatedMaxwellFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    inComingByRho_(ptf.inComingByRho_),
    outGoing_(ptf.outGoing_)
{}


template<class Type>
calculatedMaxwellFvPatchField<Type>::calculatedMaxwellFvPatchField
(
    const calculatedMaxwellFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    inComingByRho_(ptf.inComingByRho_),
    outGoing_(ptf.outGoing_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void calculatedMaxwellFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    inComingByRho_.autoMap(m);
    outGoing_.autoMap(m);
}


template<class Type>
void calculatedMaxwellFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const calculatedMaxwellFvPatchField<Type>& fgptf =
        refCast<const calculatedMaxwellFvPatchField<Type> >(ptf);

    inComingByRho_.rmap(fgptf.inComingByRho_, addr);
    outGoing_.rmap(fgptf.outGoing_, addr);
}

template<class Type>
void calculatedMaxwellFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    //- Get global outgoing density flux at wall
    //Field<Type> outGoingPart(outGoing_);
    //MPI_Allreduce(
        //outGoingPart.data(),
        //outGoing_.data(),
        //outGoing_.size(),
        //MPI_DOUBLE,
        //MPI_SUM,
        //MPI_COMM_WORLD
    //);

    fvPatchField<Type>::operator==(outGoing_/mag(inComingByRho_));
    //- Debug
    //Info << inComingByRho_ << endl;
    fvPatchField<Type>::evaluate();
    // reset outGoing to zero
    outGoing_ = pTraits<Type>::zero;
}

template<class Type>
void calculatedMaxwellFvPatchField<Type>::write(Ostream& os) const
{
        fvPatchField<Type>::write(os);
        this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
