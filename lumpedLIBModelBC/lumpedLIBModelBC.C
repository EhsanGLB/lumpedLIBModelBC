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

\*---------------------------------------------------------------------------*/

#include "lumpedLIBModelBC.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lumpedLIBModelBC::lumpedLIBModelBC
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    state_("discharge"),
    rho_(0.0),
    Cp_(0.0),
    I_(0.0),
    R_(0.0),
    dEdT_(0.0),
    TName_("T"),
    Ti_(0.0),
    SOC_(0.0),
    C_(0.0),
    V_(0.0),
    kf_(0.0)
{}


Foam::lumpedLIBModelBC::lumpedLIBModelBC
(
    const lumpedLIBModelBC& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    state_(ptf.state_),
    rho_(ptf.rho_),
    Cp_(ptf.Cp_),
    I_(ptf.I_),
    R_(ptf.R_),
    dEdT_(ptf.dEdT_),
    TName_(ptf.TName_),
    Ti_(ptf.Ti_),
    SOC_(ptf.SOC_),
    C_(ptf.C_),
    V_(ptf.V_),
    kf_(ptf.kf_)
{}


Foam::lumpedLIBModelBC::lumpedLIBModelBC
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    state_(dict.lookup("state")),
    rho_(readScalar(dict.lookup("rho"))),
    Cp_(readScalar(dict.lookup("Cp"))),
    I_(readScalar(dict.lookup("I"))),
    R_(readScalar(dict.lookup("R"))),
    dEdT_(readScalar(dict.lookup("dEdT"))),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    Ti_(readScalar(dict.lookup("Ti"))),
    SOC_(readScalar(dict.lookup("SOC"))),
    C_(readScalar(dict.lookup("C"))),
    V_(readScalar(dict.lookup("V"))),
    kf_(readScalar(dict.lookup("kf")))
{}


Foam::lumpedLIBModelBC::lumpedLIBModelBC
(
    const lumpedLIBModelBC& fcvpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(fcvpvf, iF),
    state_(fcvpvf.state_),
    rho_(fcvpvf.rho_),
    Cp_(fcvpvf.Cp_),
    I_(fcvpvf.I_),
    R_(fcvpvf.R_),
    dEdT_(fcvpvf.dEdT_),
    TName_(fcvpvf.TName_),
    Ti_(fcvpvf.Ti_),
    SOC_(fcvpvf.SOC_),
    C_(fcvpvf.C_),
    V_(fcvpvf.V_),
    kf_(fcvpvf.kf_)
{}


Foam::lumpedLIBModelBC::lumpedLIBModelBC
(
    const lumpedLIBModelBC& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    state_(ptf.state_),
    rho_(ptf.rho_),
    Cp_(ptf.Cp_),
    I_(ptf.I_),
    R_(ptf.R_),
    dEdT_(ptf.dEdT_),
    TName_(ptf.TName_),
    Ti_(ptf.Ti_),
    SOC_(ptf.SOC_),
    C_(ptf.C_),
    V_(ptf.V_),
    kf_(ptf.kf_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lumpedLIBModelBC::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const std::string red("\033[0;31m");
    const std::string reset("\033[0m");


    //- determine the sign of currernt
    scalar sign_;
    if(state_ == "discharge"){sign_ = -1.0;}
    else{sign_ = 1.0;}


    //- get the currernt time, time step, ID, and name
    scalar t_ = this->db().time().value();//this->db().time().timeOutputValue();
    scalar deltaT_ = this->db().time().deltaT().value();
    label patchID_ = this->patch().index();
    word patchName_ = this->patch().name();


    //- read names and types of patches
    const volScalarField& T_ = db().lookupObject<volScalarField>(TName_);
    const volScalarField::GeometricBoundaryField& Tbf_ = T_.boundaryField();
    wordList TBoundaryTypes_ = Tbf_.types();
    wordList patchNames_ = patch().patch().boundaryMesh().names();
    wordList patchTypes_ = patch().patch().boundaryMesh().types();
    IOList<scalar> SOCValues_( IOobject( "SOCValues", T_.mesh().time().timeName(), T_.mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE ), 20 );

    if( max(SOCValues_) <= 0.0 )
    {
        FatalErrorIn("finished") << exit(FatalError);
    }

    //- compute the SOC
    scalar SOCCurr_ = SOC_ + sign_*(100/(3600*C_))*I_*deltaT_;
    scalar Tb_ = Ti_;

    const fvPatchScalarField& Tp = lookupPatchField<volScalarField, scalar>(TName_);
    scalarField gradTp_ = Tp.snGrad();
    scalar QOhmic_ = 0.0;
    scalar QEntropic_ = 0.0;
    scalar Qb_ = 0.0;

    if( SOCCurr_ > 0.0 )
    {
        //- compute the LIB temperature by lumped model
        QOhmic_ = R_*I_*I_;
        QEntropic_ = -1.0*sign_*I_*Ti_*dEdT_;
        SOC_ = SOCCurr_;
        Info << "SOC of " << patchName_ << ": " << SOC_ << endl;
    }

    if( SOCCurr_ <= 0.0 )
    {
        //FatalErrorIn("finished") << exit(FatalError);
        QOhmic_ = 0.0;
        QEntropic_ = 0.0;
        SOC_ = 0.0;
        cout << red << "Notice: " << patchName_ << " is discharged." << reset << "\n";
    }

    Qb_ = QOhmic_ + QEntropic_;
    scalar QLoss_ = 0.0;
    forAll(Tp, i)
    {
        QLoss_ += -kf_*gradTp_[i]*patch().magSf()[i];
    }

    Tb_ = Ti_ + deltaT_*( (Qb_+QLoss_) / (rho_*Cp_*V_) );
    Info << "Temperature of " << patchName_ << ": " << Tb_ << endl;
    SOCValues_[patchID_] = SOC_;

    scalarField::operator=(Tb_);
    fixedValueFvPatchScalarField::updateCoeffs();
}


// Write
void Foam::lumpedLIBModelBC::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("state") << state_ << token::END_STATEMENT << nl;
    os.writeKeyword("rho") << rho_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cp") << Cp_ << token::END_STATEMENT << nl;
    os.writeKeyword("R") << R_ << token::END_STATEMENT << nl;
    os.writeKeyword("dEdT") << dEdT_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ti") << Ti_ << token::END_STATEMENT << nl;
    os.writeKeyword("SOC") << SOC_ << token::END_STATEMENT << nl;
    os.writeKeyword("C") << C_ << token::END_STATEMENT << nl;
    os.writeKeyword("V") << V_ << token::END_STATEMENT << nl;
    os.writeKeyword("kf") << kf_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makePatchTypeField(fvPatchScalarField, lumpedLIBModelBC);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

