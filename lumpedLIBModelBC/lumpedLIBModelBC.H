/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    lumpedLIBModelBC

Description
    Boundary condition specifies a lumped model to compute temperature of
    lithium-ion battery, given state, density, specific heat capacity, current
   , resistance, the entropy coefficient, initial temperature, state of charge
   , capacity, volume, and thermal conductivity of fluid.

    @verbatim
        inlet
        {
            type                    lumpedLIBModelBC;
            state                   discharge;
            rho                     2964;//(kg/m3)
            Cp                      1108;//(J/Kg.K)
            I                       3.2;//(A)
            R                       0.2;//(Ohm)
	    dEdT                    4.0e-4;//(V/K)
	    Ti                      293.0;//(K)
            SOC                     100;//(%)
            C                       3.2;//(Ah)
            V                       1.6532e-5;//(m3)
            kf                      0.026;//(W/m.K)
        }
    @endverbatim

SourceFiles
    lumpedLIBModelBC.C

Author
    Ehsan Golab, SUT. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef lumpedLIBModelBC_H
#define lumpedLIBModelBC_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonUniformTemperatureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class lumpedLIBModelBC
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- state
        word state_;

        //- LIB density
        scalar rho_;

        //- LIB heat specific capacity
        scalar Cp_;

        //- current
        scalar I_;

        //- resistance
        scalar R_;

        //- entropy coefficient
        scalar dEdT_;

        //- Name of the temperature field
        word TName_;

        //- initial temperature of LIB
        scalar Ti_;

        //- state of charge
        scalar SOC_;

        //- LIB capacity
        scalar C_;

        //- LIB volume
        scalar V_;

        //- thermal conductivity of fluid
        scalar kf_;


public:

    //- Runtime type information
    TypeName("lumpedLIBModelBC");


    // Constructors

        //- Construct from patch and internal field
        lumpedLIBModelBC
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        lumpedLIBModelBC
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given lumpedLIBModelBC
        //  onto a new patch
        lumpedLIBModelBC
        (
            const lumpedLIBModelBC&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new lumpedLIBModelBC(*this)
            );
        }

        //- Construct as copy setRg internal field reference
        lumpedLIBModelBC
        (
            const lumpedLIBModelBC&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setRg internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new lumpedLIBModelBC(*this, iF)
            );
        }

        //- Construct as copy
        lumpedLIBModelBC
        (
            const lumpedLIBModelBC&
        );


    // Member functions

        //- Update coefficients
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
