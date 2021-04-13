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

Foam::scalar Foam::relaxAbsorptionFvPatchVectorField::alphaR(const scalar& cellx)
{
    scalar Xi = (x0_ - cellx) / (x0_ - x1_);

    scalar alR = 1.0 - (exp(pow(Xi, 3.5)) - 1.0) / (exp(1.0) - 1.0);

    return alR;
}

void Foam::relaxAbsorptionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Initialization of zero velocity field at the boundary
    vectorField boundaryU(this->patch().size(), vector::zero);

    // forAll (boundaryU, facei)
    // {
    //     boundaryU[facei] = vector(uInterp, 0, wInterp);
    // }

    // Apply relaxation zone for the internal mesh
    const volVectorField& U = db().lookupObject<volVectorField>("U");
    volVectorField& Ucast = const_cast<volVectorField&>(U);

    forAll(U.mesh().cells(),celli)
    {
        const scalar& cellx = U.mesh().C()[celli][0];

        if ((cellx >= x0_) && (cellx <= x1_))
        {
            Ucast.primitiveFieldRef()[celli] = U[celli] * alphaR(cellx);
        }
        else if (cellx > x1_)
        {
            Ucast.primitiveFieldRef()[celli] = vector(0.0, 0.0, 0.0);
        }
    }

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
