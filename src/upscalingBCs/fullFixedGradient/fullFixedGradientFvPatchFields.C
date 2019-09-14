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

#include "fullFixedGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"



// for non-templated patch fields
#define mymakePatchTypeField(PatchTypeField, typePatchTypeField)                \
    defineTypeNameAndDebug(typePatchTypeField, 0);                            \
    addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)

// for non-templated patch fields - use with caution
#define mymakeRemovablePatchTypeField(PatchTypeField, typePatchTypeField)       \
    defineTypeNameAndDebug(typePatchTypeField, 0);                            \
    addRemovableToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)


// for templated patch fields
#define mymakeTemplatePatchTypeField(PatchTypeField, typePatchTypeField)        \
    defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);               \
    addToPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)


#define mymakePatchFields(type)                                                 \
    mymakeTemplatePatchTypeField                                                \
    (                                                                         \
        fvPatchScalarField,                                                   \
        type##FvPatchScalarField                                              \
    );                                                                        \
    mymakeTemplatePatchTypeField                                                \
    (                                                                         \
        fvPatchVectorField,                                                   \
        type##FvPatchVectorField                                              \
    );                                                                        \


#define mymakePatchFieldsTypeName(type)                                         \
    defineNamedTemplateTypeNameAndDebug(type##FvPatchScalarField, 0);         \
    defineNamedTemplateTypeNameAndDebug(type##FvPatchVectorField, 0);         \


#define mymakePatchTypeFieldTypedefs(type)                                      \
    typedef type##FvPatchField<scalar> type##FvPatchScalarField;              \
    typedef type##FvPatchField<vector> type##FvPatchVectorField;              \



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

mymakePatchFields(fullFixedGradient);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
