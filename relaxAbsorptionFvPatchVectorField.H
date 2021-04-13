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

Class
    Foam::relaxAbsorptionFvPatchVectorField

Description
    Wave absorption using relaxation zone

SourceFiles
    relaxAbsorptionFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef relaxAbsorptionFvPatchVectorField_H
#define relaxAbsorptionFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
         Class relaxAbsorptionFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class relaxAbsorptionFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Begin of the relaxation zone
        scalar x0_;

        //- End of the relaxation zone
        scalar x1_;


public:

   //- Runtime type information
   TypeName("relaxAbsorption");


   // Constructors

        //- Construct from patch and internal field
        relaxAbsorptionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        relaxAbsorptionFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given relaxAbsorptionFvPatchVectorField
        //  onto a new patch
        relaxAbsorptionFvPatchVectorField
        (
            const relaxAbsorptionFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );


    // Member functions
	
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
