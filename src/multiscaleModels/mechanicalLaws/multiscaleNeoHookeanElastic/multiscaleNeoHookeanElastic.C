/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "multiscaleNeoHookeanElastic.H"
#include "addToRunTimeSelectionTable.H"
#include <thread>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiscaleNeoHookeanElastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, multiscaleNeoHookeanElastic, nonLinGeomMechLaw
    );
}

// * * * * * * * * * * *  Private Member Funtcions * * * * * * * * * * * * * //

void Foam::multiscaleNeoHookeanElastic::makeF()
{
    if (FPtr_)
    {
        FatalErrorIn("void Foam::multiscaleNeoHookeanElastic::makeF()")
            << "pointer already set" << abort(FatalError);
    }

    FPtr_ =
        new volTensorField
        (
            IOobject
            (
                "F",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );
}


Foam::volTensorField& Foam::multiscaleNeoHookeanElastic::F()
{
    if (!FPtr_)
    {
        makeF();
    }

    return *FPtr_;
}


void Foam::multiscaleNeoHookeanElastic::makeFf()
{
    if (FfPtr_)
    {
        FatalErrorIn("void Foam::multiscaleNeoHookeanElastic::makeFf()")
            << "pointer already set" << abort(FatalError);
    }

    FfPtr_ =
        new surfaceTensorField
        (
            IOobject
            (
                "Ff",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );
}


Foam::surfaceTensorField& Foam::multiscaleNeoHookeanElastic::Ff()
{
    if (!FfPtr_)
    {
        makeFf();
    }

    return *FfPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::multiscaleNeoHookeanElastic::multiscaleNeoHookeanElastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    applyBoundaryRve_
    (
         dict.lookupOrDefault<bool>("applyBoundaryRve", true)
    ),
    pRveSolvers_(),
    rho_(dict.lookup("rho")),
    E_(dict.lookup("E")),
    nu_(dict.lookup("nu")),
    mu_(E_/(2.0*(1.0 + nu_))),
    K_
    (
        planeStress()
      ? (nu_*E_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*mu_
      : (nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*mu_
    ),
    pE_(E_*0.4),
    pmu_(pE_/(2.0*(1.0 + nu_))),
    pK_
    (
        planeStress()
      ? (nu_*pE_/((1.0 + nu_)*(1.0 - nu_))) + (2.0/3.0)*pmu_
      : (nu_*pE_/((1.0 + nu_)*(1.0 - 2.0*nu_))) + (2.0/3.0)*pmu_
    ),
    FPtr_(NULL),
    FfPtr_(NULL)
{
    F().oldTime();
  
    // Determine number of rve (cells + non-empty boundary faces)
    label nRve = mesh.nCells();

    if (applyBoundaryRve_)
    {
    forAll(mesh.boundary(), patchI)
    {
        if (!mesh.boundary()[patchI].coupled())
        {
            nRve += mesh.boundary()[patchI].size();
        }
    }
    }

    Pout << "Number of RVEs : " << nRve << endl;

    // Create periodic RVE solvers
    pRveSolvers_.setSize(nRve);
    forAll(pRveSolvers_, nRve)
    {
        OStringStream RegionName;
        RegionName() << "rve-" << nRve;

        if (Pstream::parRun())
        {
            fileName src = mesh.time().path()/".."/".."/"rve-periodic"/"constant";
            fileName dst = mesh.time().path()/"constant"/word(RegionName.str());
            ln(src, dst);
        
            //src = mesh.time().path()/".."/".."/"rve-periodic"/"constant"/"plasticStrainVsYieldStress";
            //dst = mesh.time().path()/".."/"constant"/"plasticStrainVsYieldStress";
            //ln(src, dst);
            
            src = mesh.time().path()/".."/".."/"rve-periodic"/"constant";
            dst = mesh.time().path()/".."/"constant"/word(RegionName.str());
            ln(src, dst);
            
            src = mesh.time().path()/".."/".."/"rve-periodic"/"0";
            dst = mesh.time().path()/"0"/word(RegionName.str());
            ln(src, dst);
        
            src = mesh.time().path()/".."/".."/"rve-periodic"/"system";
            dst = mesh.time().path()/".."/"system"/word(RegionName.str());
            ln(src, dst);
        }
        else
        {
            fileName src = mesh.time().path()/".."/"rve-periodic"/"constant";
            fileName dst = mesh.time().path()/"constant"/word(RegionName.str());
            ln(src, dst);
        
            src = mesh.time().path()/".."/"rve-periodic"/"constant"/"plasticStrainVsYieldStress";
            dst = mesh.time().path()/"constant"/"plasticStrainVsYieldStress";
            ln(src, dst);
            
            src = mesh.time().path()/".."/"rve-periodic"/"0";
            dst = mesh.time().path()/"0"/word(RegionName.str());
            ln(src, dst);
        
            src = mesh.time().path()/".."/"rve-periodic"/"system";
            dst = mesh.time().path()/"system"/word(RegionName.str());
            ln(src, dst);
        }

        pRveSolvers_.set
        (
            nRve,

            new solidModels::pRveNonLinGeomTotalLagSolid

            (
                const_cast<Time&>(mesh.time()),
                word(RegionName.str())
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiscaleNeoHookeanElastic::~multiscaleNeoHookeanElastic()
{
    deleteDemandDrivenData(FPtr_);
    deleteDemandDrivenData(FfPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::multiscaleNeoHookeanElastic::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "rhoLaw",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            rho_,
            calculatedFvPatchScalarField::typeName
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::multiscaleNeoHookeanElastic::impK() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "impK",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            (4.0/3.0)*mu_ + K_ // == 2*mu + lambda
        )
    );
}


void Foam::multiscaleNeoHookeanElastic::correct(volSymmTensorField& sigma)
{
    // Check if the mathematical model is in total or updated Lagrangian form
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Calculate the relative deformation gradient
        const volTensorField relF = I + gradDD.T();

        // Update the total deformation gradient
        F() = relF & F().oldTime();
    }
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        // Lookup gradient of displacement
        const volTensorField& gradD =
            mesh().lookupObject<volTensorField>("grad(D)");

        // Update the total deformation gradient
        F() = I + gradD.T();

        // Update the relative deformation gradient: not needed
        volTensorField relF = F() & inv(F().oldTime());

        // Set relative deformation gradient for pRve solvers
        forAll(F().internalField(), cellI)
        {
            
                pRveSolvers_[cellI].avgRelF() = relF.internalField()[cellI];
        }
        
        if (applyBoundaryRve_)
        {
        //label solI = mesh().nCells();
        //label nRve = mesh().nCells();
        //Pout << "Number of RVEs 1 : " << nRve << endl;
        //forAll(relF.boundaryField(), patchI)
        //{
          //  forAll(relF.boundaryField()[patchI],faceI)
            //{
              //      pRveSolvers_[nRve].avgRelF() =
                //      relF.boundaryField()[patchI][faceI];
            //nRve++;
            //}
        //}
        }
        //Pout << "Number of RVEs 2 : " << nRve << endl;

        // Solve pRve
        forAll(pRveSolvers_, nRve)
        {
            Info << "Solving pRve: " << nRve << endl;
            pRveSolvers_[nRve].evolve();
        }
        
        // Get average stess from RVE model
        forAll(sigma.internalField(), cellI)
        {
                    sigma.internalField()[cellI] =
                    pRveSolvers_[cellI].avgSigma().value();
        }
        
        if (applyBoundaryRve_)
        {
        //solI = mesh().nCells();
        //label nRve = mesh().nCells();
        //forAll(sigma.boundaryField(), patchI)
        //{
          //      forAll(sigma.boundaryField()[patchI], faceI)
            //    {
              //      sigma.boundaryField()[patchI][faceI] =
                //        pRveSolvers_[nRve].avgSigma().value();
                  //  nRve++;
                    //solI++;
                //}
        //}
        }
        else
        {
            // Calculate the Jacobian of the deformation gradient
            forAll(sigma.boundaryField(), patchI)
            {
                const scalarField pJ =
                    det(F().boundaryField()[patchI]);

                const tensorField pF =
                    F().boundaryField()[patchI];

                //Calculate the volume preserving left Cauchy Green strain
                const symmTensorField pbEbar =
                    pow(pJ, -2.0/3.0)*symm(pF & pF.T());

                //Calculate the deviatoric stress
                const symmTensorField ps = mu_.value()*dev(pbEbar);

                //Calculate the Cauchy stress
                sigma.boundaryField()[patchI] =
                    (1.0/pJ)*(0.5*K_.value()*(pow(pJ, 2) - 1)*I + ps);
            }
        }

        sigma.correctBoundaryConditions();
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::multiscaleNeoHookeanElastic::"
            "correct(volSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }
}


void Foam::multiscaleNeoHookeanElastic::correct(surfaceSymmTensorField& sigma)
{
    // Check if the mathematical model is in total or updated Lagrangian form
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        // Update the relative deformation gradient: not needed
        const surfaceTensorField relF = I + gradDD.T();

        // Update the total deformation gradient
        Ff() = relF & Ff().oldTime();
    }
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
        // Lookup gradient of displacement
        const surfaceTensorField& gradD =
            mesh().lookupObject<surfaceTensorField>("grad(D)f");

        // Update the total deformation gradient
        Ff() = I + gradD.T();

        // Update the relative deformation gradient: not needed
        //relF() = F() & inv(F().oldTime());
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::multiscaleNeoHookeanElastic::"
            "correct(surfaceSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }
}

void Foam::multiscaleNeoHookeanElastic::updateTotalFields()
{
    forAll(pRveSolvers_, nRve)
    {
        pRveSolvers_[nRve].updateTotalFields();
    }
}

// ************************************************************************* //
