/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "RobinFvPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::RobinFvPatchField<Type>::RobinFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    RobinD_(p.size()),
    RobinK_(p.size()),
    RobinF_(p.size()),
    dt_("DT"),
    drift_("none"),
    FField_("none")
{
}


template<class Type>
Foam::RobinFvPatchField<Type>::RobinFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict),
    RobinD_("RobinD", dict, p.size()),
    RobinK_("RobinK", dict, p.size()),
    RobinF_("RobinF", dict, p.size()),
    dt_(dict.lookupOrDefault<word>("dt","none")),
    drift_(dict.lookupOrDefault<word>("drift","none")),
    FField_(dict.lookupOrDefault<word>("FField","none"))
{
    evaluate();
}


template<class Type>
Foam::RobinFvPatchField<Type>::RobinFvPatchField
(
    const RobinFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
    RobinD_(ptf.RobinD_, mapper),
    RobinK_(ptf.RobinK_, mapper),
    RobinF_(ptf.RobinF_, mapper),
    dt_(ptf.dt_),
    drift_(ptf.drift_),
    FField_(ptf.FField_)
{
    if (notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::RobinFvPatchField<Type>::RobinFvPatchField
(
    const RobinFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    RobinD_(ptf.RobinD_),
    RobinK_(ptf.RobinK_),
    RobinF_(ptf.RobinF_),
    dt_(ptf.dt_),
    drift_(ptf.drift_),
    FField_(ptf.FField_)
{}


template<class Type>
Foam::RobinFvPatchField<Type>::RobinFvPatchField
(
    const RobinFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    RobinD_(ptf.RobinD_),
    RobinK_(ptf.RobinK_),
    RobinF_(ptf.RobinF_),
    dt_(ptf.dt_),
    drift_(ptf.drift_),
    FField_(ptf.FField_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::RobinFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    RobinD_.autoMap(m);
    RobinK_.autoMap(m);
    RobinF_.autoMap(m);
    //drift_.autoMap(m);
}


template<class Type>
void Foam::RobinFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const RobinFvPatchField<Type>& mptf =
        refCast<const RobinFvPatchField<Type>>(ptf);

    RobinD_.rmap(mptf.RobinD_, addr);
    RobinK_.rmap(mptf.RobinK_, addr);
    RobinF_.rmap(mptf.RobinF_, addr);
    //drift_.rmap(mptf.drift_, addr);
}


template<class Type>
void Foam::RobinFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        (
          this->patchInternalField() * RobinD2_ * this->patch().deltaCoeffs()
          + RobinF2_
        )
        /
        (
           RobinD2_ * this->patch().deltaCoeffs()
          - RobinK2_
        )
    );

    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::snGrad() const
{
    return
        (
          RobinK2_ * this->patchInternalField()
          + RobinF2_
        ) * this->patch().deltaCoeffs()
        /
        (
          RobinD2_ * this->patch().deltaCoeffs()
         - RobinK2_
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return Type(pTraits<Type>::one) *
    ( 1.0 + RobinK2_/(RobinD2_*this->patch().deltaCoeffs() - RobinK2_) );

}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
         (
           RobinF2_         /
           (
            RobinD2_ * this->patch().deltaCoeffs()
           - RobinK2_
           )
         );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::gradientInternalCoeffs() const
{
    return Type(pTraits<Type>::one) *
    (
      RobinK2_*this->patch().deltaCoeffs()
      /
      (
        RobinD2_ * this->patch().deltaCoeffs()
       - RobinK2_
      )
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::RobinFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
       (
         RobinF2_ * this->patch().deltaCoeffs()
         /
         (
          RobinD2_ * this->patch().deltaCoeffs()
         - RobinK2_
         )
       );
}


template<class Type>
void Foam::RobinFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    RobinD_.writeEntry("RobinD", os);
    RobinK_.writeEntry("RobinK", os);
    RobinF_.writeEntry("RobinF", os);
    os.writeKeyword("drift") << drift_ << token::END_STATEMENT << nl;
    os.writeKeyword("FField") << FField_ << token::END_STATEMENT << nl;
    os.writeKeyword("dt") << dt_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}



template<class Type>
void Foam::RobinFvPatchField<Type>::updateCoeffs()
 {
     if (this->updated())
     {
         return;
     }

     //const IOdictionary& transportProperties = this->db().template lookupObject<IOdictionary>("transportProperties");
     //const dimensionedScalar DT(transportProperties.lookup("DT"));
//        const GeometricField<Type, fvsPatchField, surfaceMesh> gradT(fvc::snGrad(
//            this->db().template
//            lookupObject<GeometricField<Type, fvPatchField, volMesh> >(fieldName_)));
     if (dt_ != "none")
     {
        const IOdictionary& transportProperties = this->db().template lookupObject<IOdictionary>("transportProperties");
        const dimensionedTensor DT(transportProperties.lookup("DT"));
        const tensorField DTF(RobinD_.size(), DT.value());
        //const dimensionedTensor& DT
       // (
       //  this->db().template lookupObject<dimensionedScalar>(dt_)
       // );
      RobinD2_ = ((DTF & this->patch().nf()) & this->patch().nf()) + RobinD_;
     }
     else
     {
       RobinD2_ = RobinD_;
     }

     if (drift_ != "none")
     {
       const fvPatchVectorField& Up
        (
         this->patch().template lookupPatchField<volVectorField, vector>(drift_)
        );
        RobinK2_ = RobinK_ + (Up & this->patch().nf());
     }
     else
     {
       RobinK2_ = RobinK_;
     }

     if (FField_ != "none")
     {
       const fvPatchField<Type>& FF
        (
         this->patch().template lookupPatchField<GeometricField<Type, fvPatchField, volMesh>, Type>(FField_)
        );
        RobinF2_ = RobinF_ + FF;
     }
     else
     {
       RobinF2_ = RobinF_;
     }

     fvPatchField<Type>::updateCoeffs();
 }

// ************************************************************************* //
