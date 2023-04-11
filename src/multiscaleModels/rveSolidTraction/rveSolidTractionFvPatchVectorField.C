/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "rveSolidTractionFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "lookupSolidModel.H"
#include "pRveNonLinGeomTotalLagSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rveSolidTractionFvPatchVectorField::
rveSolidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_(),
    secondOrder_(false),
    limitCoeff_(1.0),
    relaxFac_(1.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = vector::zero;
}


rveSolidTractionFvPatchVectorField::
rveSolidTractionFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    traction_(p.size(), vector::zero),
    pressure_(p.size(), 0.0),
    tractionSeries_(),
    pressureSeries_(),
    secondOrder_(dict.lookupOrDefault<Switch>("secondOrder", false)),
    limitCoeff_(dict.lookupOrDefault<scalar>("limitCoeff", 1.0)),
    relaxFac_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0))
{
    Info<< "Creating " << type() << " boundary condition" << endl;

    if (dict.found("gradient"))
    {
        gradient() = vectorField("gradient", dict, p.size());
    }
    else
    {
        gradient() = vector::zero;
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        fvPatchVectorField::operator=(patchInternalField());
    }

    // Check if traction is time-varying
    if (dict.found("tractionSeries"))
    {
        Info<< "    traction is time-varying" << endl;
        tractionSeries_ =
            interpolationTable<vector>(dict.subDict("tractionSeries"));
    }
    else
    {
        traction_ = vectorField("traction", dict, p.size());
    }

    // Check if pressure is time-varying
    if (dict.found("pressureSeries"))
    {
        Info<< "    pressure is time-varying" << endl;
        pressureSeries_ =
            interpolationTable<scalar>(dict.subDict("pressureSeries"));
    }
    else
    {
        pressure_ = scalarField("pressure", dict, p.size());
    }

    if (secondOrder_)
    {
        Info<< "    second order correction" << endl;
    }

    if (limitCoeff_)
    {
        Info<< "    limiter coefficient: " << limitCoeff_ << endl;
    }

    if (relaxFac_ < 1.0)
    {
        Info<< "    relaxation factor: " << relaxFac_ << endl;
    }
}


rveSolidTractionFvPatchVectorField::
rveSolidTractionFvPatchVectorField
(
    const rveSolidTractionFvPatchVectorField& stpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(stpvf, p, iF, mapper),
    traction_(stpvf.traction_, mapper),
    pressure_(stpvf.pressure_, mapper),
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_),
    secondOrder_(stpvf.secondOrder_),
    limitCoeff_(stpvf.limitCoeff_),
    relaxFac_(stpvf.relaxFac_)
{}


rveSolidTractionFvPatchVectorField::
rveSolidTractionFvPatchVectorField
(
    const rveSolidTractionFvPatchVectorField& stpvf
)
:
    fixedGradientFvPatchVectorField(stpvf),
    traction_(stpvf.traction_),
    pressure_(stpvf.pressure_),
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_),
    secondOrder_(stpvf.secondOrder_),
    limitCoeff_(stpvf.limitCoeff_),
    relaxFac_(stpvf.relaxFac_)
{}


rveSolidTractionFvPatchVectorField::
rveSolidTractionFvPatchVectorField
(
    const rveSolidTractionFvPatchVectorField& stpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(stpvf, iF),
    traction_(stpvf.traction_),
    pressure_(stpvf.pressure_),
    tractionSeries_(stpvf.tractionSeries_),
    pressureSeries_(stpvf.pressureSeries_),
    secondOrder_(stpvf.secondOrder_),
    limitCoeff_(stpvf.limitCoeff_),
    relaxFac_(stpvf.relaxFac_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void rveSolidTractionFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    traction_.autoMap(m);
    pressure_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void rveSolidTractionFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const rveSolidTractionFvPatchVectorField& dmptf =
        refCast<const rveSolidTractionFvPatchVectorField>(ptf);

    traction_.rmap(dmptf.traction_, addr);
    pressure_.rmap(dmptf.pressure_, addr);
}


// Update the coefficients associated with the patch field
void rveSolidTractionFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (tractionSeries_.size())
    {
        traction_ = tractionSeries_(this->db().time().timeOutputValue());
    }

    if (pressureSeries_.size())
    {
        pressure_ = pressureSeries_(this->db().time().timeOutputValue());
    }

    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(patch().boundaryMesh().mesh());

    if (isA<solidModels::pRveNonLinGeomTotalLagSolid >(solMod))
    {
        solidModels::pRveNonLinGeomTotalLagSolid& pRveSolMod =
            const_cast<solidModels::pRveNonLinGeomTotalLagSolid&>
            (
                refCast<const solidModels::pRveNonLinGeomTotalLagSolid>(solMod)
            );

        pRveSolMod.updateAvgSigma();
        
        const dimensionedSymmTensor& curAvgSigma = pRveSolMod.avgSigma();

        const vectorField nf = patch().nf();

        traction_ = (nf & curAvgSigma.value());
        pressure_ = 0;
    }
    else
    {
        FatalErrorIn
        (
            "void rveSolidTractionFvPatchVectorField::updateCoeffs()"
        )   << "Please use solid model: "
             << solidModels::pRveNonLinGeomTotalLagSolid::typeName
             << " with this boundary condition."
             << abort(FatalError);
    }
      
    
    // Set surface-normal gradient on the patch corresponding to the desired
    // traction
    gradient() =
        relaxFac_*solMod.tractionBoundarySnGrad
        (
            traction_, pressure_, patch()
        )
      + (1.0 - relaxFac_)*gradient();

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void rveSolidTractionFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    // Lookup the gradient field
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + dimensionedInternalField().name() + ")"
        );

    // Face unit normals
    const vectorField n = patch().nf();

    // Delta vectors
    const vectorField delta = patch().delta();

    // Non-orthogonal correction vectors
    const vectorField k = ((I - sqr(n)) & delta);

    if (secondOrder_)
    {
        const vectorField dUP = (k & gradField.patchInternalField());
        const vectorField nGradUP = (n & gradField.patchInternalField());

        Field<vector>::operator=
        (
            patchInternalField()
          + dUP
          + 0.5*(gradient() + nGradUP)/patch().deltaCoeffs()
        );
    }
    else
    {

        Field<vector>::operator=
        (
            patchInternalField()
          + (k & gradField.patchInternalField())
          + gradient()/patch().deltaCoeffs()
        );
    }

    fvPatchField<vector>::evaluate();
}


void rveSolidTractionFvPatchVectorField::write(Ostream& os) const
{
    fixedGradientFvPatchVectorField::write(os);

    if (tractionSeries_.size())
    {
        os.writeKeyword("tractionSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        tractionSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        traction_.writeEntry("traction", os);
    }

    if (pressureSeries_.size())
    {
        os.writeKeyword("pressureSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
        pressureSeries_.write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        pressure_.writeEntry("pressure", os);
    }
    os.writeKeyword("secondOrder")
        << secondOrder_ << token::END_STATEMENT << nl;
    os.writeKeyword("limitCoeff")
        << limitCoeff_ << token::END_STATEMENT << nl;
    os.writeKeyword("relaxationFactor")
        << relaxFac_ << token::END_STATEMENT << nl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, rveSolidTractionFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
