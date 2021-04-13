/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenFOAM Foundation
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

#include "relaxAbsorptionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::relaxAbsorptionFvPatchVectorField::relaxAbsorptionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    x0_(0.0),
    x1_(0.0)
{}


Foam::relaxAbsorptionFvPatchVectorField::relaxAbsorptionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    x0_(readScalar(dict.lookup("x0"))),
    x1_(readScalar(dict.lookup("x1")))
{}


Foam::relaxAbsorptionFvPatchVectorField::relaxAbsorptionFvPatchVectorField
(
    const relaxAbsorptionFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    x0_(ptf.x0_),
    x1_(ptf.x1_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::relaxAbsorptionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Alternative initialization of the velocity field
    vectorField boundaryU(this->patch().size(), vector::zero);

    // forAll (boundaryU, facei)
    // {
    //     boundaryU[facei] = vector(uInterp, 0, wInterp);
    // }

    // Set zero velocity boundary field
    operator==(boundaryU);

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::relaxAbsorptionFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry("x0", x0_);
    os.writeEntry("x1", x1_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       relaxAbsorptionFvPatchVectorField
   );
}


// ************************************************************************* //
