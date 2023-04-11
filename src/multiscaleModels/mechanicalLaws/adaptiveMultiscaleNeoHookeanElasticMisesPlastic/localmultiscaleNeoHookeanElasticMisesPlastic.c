/*---------------------------------------------------------------------------*\
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

#include "localmultiscaleNeoHookeanElasticMisesPlastic.H"
#include "addToRunTimeSelectionTable.H"
#include "transformGeometricField.H"
#include "logVolFields.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(localmultiscaleNeoHookeanElasticMisesPlastic, 0);
    addToRunTimeSelectionTable
    (
        mechanicalLaw, localmultiscaleNeoHookeanElasticMisesPlastic, nonLinGeomMechLaw
    );

// * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * * //

    // Tolerance for Newton loop
    scalar localmultiscaleNeoHookeanElasticMisesPlastic::LoopTol_ = 1e-8;

    // Maximum number of iterations for Newton loop
    label localmultiscaleNeoHookeanElasticMisesPlastic::MaxNewtonIter_ = 200;

    // finiteDiff is the delta for finite difference differentiation
    scalar localmultiscaleNeoHookeanElasticMisesPlastic::finiteDiff_ = 0.25e-6;

    // Store sqrt(2/3) as we use it often
    scalar localmultiscaleNeoHookeanElasticMisesPlastic::sqrtTwoOverThree_ = ::sqrt(2.0/3.0);

} // End of namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::makeRelF()
{
    if (relFPtr_)
    {
        FatalErrorIn("void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::makeRelF()")
            << "pointer already set" << abort(FatalError);
    }

    relFPtr_ =
        new volTensorField
        (
            IOobject
            (
                "lawRelF",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );
}


Foam::volTensorField& Foam::localmultiscaleNeoHookeanElasticMisesPlastic::relF()
{
    if (!relFPtr_)
    {
        makeRelF();
    }

    return *relFPtr_;
}


void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::makeRelFf()
{
    if (relFfPtr_)
    {
        FatalErrorIn("void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::makeRelFf()")
            << "pointer already set" << abort(FatalError);
    }

    relFfPtr_ =
        new surfaceTensorField
        (
            IOobject
            (
                "lawRelFf",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );
}


Foam::surfaceTensorField& Foam::localmultiscaleNeoHookeanElasticMisesPlastic::relFf()
{
    if (!relFfPtr_)
    {
        makeRelFf();
    }

    return *relFfPtr_;
}


void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::makeJ()
{
    if (JPtr_)
    {
        FatalErrorIn("void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::makeJ()")
            << "pointer already set" << abort(FatalError);
    }

    JPtr_ =
        new volScalarField
        (
            IOobject
            (
                "lawJ",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        );

    // Store the old-time
    JPtr_->oldTime();
}


Foam::volScalarField& Foam::localmultiscaleNeoHookeanElasticMisesPlastic::J()
{
    if (!JPtr_)
    {
        makeJ();
    }

    return *JPtr_;
}


void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::makeJf()
{
    if (JfPtr_)
    {
        FatalErrorIn("void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::makeJf()")
            << "pointer already set" << abort(FatalError);
    }

    JfPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "lawJf",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1.0)
        );

    // Store the old-time
    JfPtr_->oldTime();
}


Foam::surfaceScalarField& Foam::localmultiscaleNeoHookeanElasticMisesPlastic::Jf()
{
    if (!JfPtr_)
    {
        makeJf();
    }

    return *JfPtr_;
}


void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::makeF()
{
    if (FPtr_)
    {
        FatalErrorIn("void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::makeF()")
            << "pointer already set" << abort(FatalError);
    }

    FPtr_ =
        new volTensorField
        (
            IOobject
            (
                "lawF",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );

    // Store the old-time
    FPtr_->oldTime();
}


Foam::volTensorField& Foam::localmultiscaleNeoHookeanElasticMisesPlastic::F()
{
    if (!FPtr_)
    {
        makeF();
    }

    return *FPtr_;
}


void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::makeFf()
{
    if (FfPtr_)
    {
        FatalErrorIn("void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::makeFf()")
            << "pointer already set" << abort(FatalError);
    }

    FfPtr_ =
        new surfaceTensorField
        (
            IOobject
            (
                "lawFf",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh(),
            dimensionedTensor("I", dimless, I)
        );

    // Store the old-time
    FfPtr_->oldTime();
}


Foam::surfaceTensorField& Foam::localmultiscaleNeoHookeanElasticMisesPlastic::Ff()
{
    if (!FfPtr_)
    {
        makeFf();
    }

    return *FfPtr_;
}


Foam::scalar Foam::localmultiscaleNeoHookeanElasticMisesPlastic::curYieldStress
(
    const scalar curEpsilonPEq,    // Current equivalent plastic strain
    const scalar J                 // Current Jacobian
) const
{
    // We assume that the stress-strain curve was specifed as Cauchy stress vs
    // true strain, but we want the Kirchhoff (tau) yield stress,
    // so we multiply Cauchy stress by J as tauSigmaY = J*sigmaCauchySigmaY

    return J*stressPlasticStrainSeries_(max(curEpsilonPEq, SMALL));
}


Foam::scalar Foam::localmultiscaleNeoHookeanElasticMisesPlastic::yieldFunction
(
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar DLambda,          // Plastic multiplier
    const scalar muBar,            // Scaled shear modulus
    const scalar J                 // Current Jacobian
) const
{
    // Evaluate current yield function
    // fy = mag(s) - sqrt(2/3)*curSigmaY
    // fy = mag(sTrial - 2*muBar*DLambda*plasticN) - ::sqrt(2.0/3.0)*curSigmaY;
    // fy = magSTrial - 2*muBar*DLambda - ::sqrt(2.0/3.0)*curSigmaY;
    // where
    // fy is the current value of the yield function - zero at convergence.
    // s is the current deviatoric component of tau
    // sTrial is trial version of s
    // plasticN is the return direction
    // DLambda is the current increment of plastic strain multiplier
    // curSigmaY is the current Kirchhoff yield stress which is typically a
    // function of total equivalent plastic strain (epsilonPEq + DEpsilonPEq)

    return
        magSTrial - 2*muBar*DLambda
      - sqrtTwoOverThree_
           *curYieldStress
            (
                epsilonPEqOld + sqrtTwoOverThree_*DLambda,
                J
            );
}


void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::newtonLoop
(
    scalar& DLambda,               // Plastic multiplier
    scalar& curSigmaY,             // Current yield stress
    const scalar epsilonPEqOld,    // Old equivalent plastic strain
    const scalar magSTrial,        // Deviatoric trial stress magnitude
    const scalar muBar,            // Scaled shear modulus
    const scalar J,                // Current Jacobian
    const scalar maxMagDEpsilon    // Max strain increment magnitude
) const
{
    // Loop to determine DEpsilonPEq
    // using Newton's method

    int i = 0;
    scalar fTrial = yieldFunction(epsilonPEqOld, magSTrial, DLambda, muBar, J);
    scalar residual = 1.0;
    do
    {
        // Numerically calculate derivative of yield function fTrial

        // First Order
        // Hauser 2009 suggested First Order is more efficient for Newton's
        // Method as only two function evaluations are required.

        // fTrial step is the the value of fTrial after a small finite
        // difference step
        const scalar fTrialStep  =
            yieldFunction
            (
                epsilonPEqOld, magSTrial, DLambda + finiteDiff_, muBar, J
            );

        // Numerical derivative of fTrial
        const scalar fTrialDerivative = (fTrialStep - fTrial)/finiteDiff_;

        // Update DLambda
        residual = fTrial/fTrialDerivative;
        DLambda -= residual;

        residual /= maxMagDEpsilon; // Normalise wrt max strain increment

        // fTrial will go to zero at convergence
        fTrial = yieldFunction(epsilonPEqOld, magSTrial, DLambda, muBar,  J);

        if (i == MaxNewtonIter_)
        {
            WarningIn("localmultiscaleNeoHookeanElasticMisesPlastic::newtonLoop()")
                << "Plasticity Newton loop not converging" << endl;
        }
    }
    while ((mag(residual) > LoopTol_) && ++i < MaxNewtonIter_);

    // Update current yield stress
    // Note: we divide by J to change the Kirchhoff yield stress to Cauchy yield
    // stress
    curSigmaY =
        curYieldStress
        (
            epsilonPEqOld + sqrtTwoOverThree_*DLambda,  J
        )/J;
}


Foam::tmp<Foam::volScalarField> Foam::localmultiscaleNeoHookeanElasticMisesPlastic::Ibar
(
    const volSymmTensorField& devBEbar
)
{
    // From Simo & Hughes 1998:
    // but this incorrectly results in det(bEbar) =! 1
    //bEbar = (s/mu) + Ibar*I;

    // A method of calculating Ibar to enforce det(bEbar) == 1 is proposed
    // by solving a cubic equation.
    // Rubin and Attia, CALCULATION OF HYPERELASTIC RESPONSE OF FINITELY
    // DEFORMED ELASTIC-VISCOPLASTIC MATERIALS, INTERNATIONAL JOURNAL FOR
    // NUMERICAL METHODS IN ENGINEERING, VOL. 39,309-320(1996)
    // and
    // M. Hollenstein M. Jabareen M. B. Rubin, Modeling a smooth elastic-
    // inelastic transition with a strongly objective numerical integrator
    // needing no iteration, Comput Mech (2013) 52:649–667
    // DOI 10.1007/s00466-013-0838-7

    // Note: In Hollenstrain et al. (2013), they suggest that Eqn 49(a) in the
    // original Rubin and Attia paper should be used.

    // Method implemented below is translated from the SmoothMultiPhase fortran
    // subroutine of Rubin

    tmp<volScalarField> tIbar
    (
        new volScalarField
        (
            IOobject
            (
                "Ibar",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    volScalarField& Ibar = tIbar();

    // Take reference to internal fields for efficiency
    scalarField& IbarI = Ibar.internalField();
    const symmTensorField devBEbarI = devBEbar.internalField();

    // Calculate internal field
    forAll(IbarI, cellI)
    {
        const scalar detdevBepr = det(devBEbarI[cellI]);
        const scalar dotprod = devBEbarI[cellI] && devBEbarI[cellI];
        const scalar fac1 = 2.0*dotprod/3.0;

        scalar alpha1 = 0.0;

        if (mag(fac1) < SMALL)
        {
            alpha1 = 3.0;
        }
        else
        {
            const scalar fac2 = (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

            if (fac2 >= 1.0)
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cosh(Foam::acosh(fac2)/3.0);
            }
            else
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cos(Foam::acos(fac2)/3.0);
            }
        }

        IbarI[cellI] = alpha1/3.0;
    }

    // Calculate boundary field
    forAll(Ibar.boundaryField(), patchI)
    {
        if
        (
            !Ibar.boundaryField()[patchI].coupled()
         && Ibar.boundaryField()[patchI].type() != "empty"
        )
        {
            // Take reference to patch fields for efficiency
            scalarField& IbarP = Ibar.boundaryField()[patchI];
            const symmTensorField& devBEbarP =
                devBEbar.boundaryField()[patchI];

            forAll(IbarP, faceI)
            {
                const scalar detdevBepr = det(devBEbarP[faceI]);
                const scalar dotprod =
                    devBEbarP[faceI] && devBEbarP[faceI];
                const scalar fac1 = 2.0*dotprod/3.0;

                scalar alpha1 = 0.0;

                if (mag(fac1) < SMALL)
                {
                    alpha1 = 3.0;
                }
                else
                {
                    const scalar fac2 =
                        (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

                    if (fac2 >= 1.0)
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cosh(Foam::acosh(fac2)/3.0);
                    }
                    else
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cos(Foam::acos(fac2)/3.0);
                    }
                }

                IbarP[faceI] = alpha1/3.0;
            }
        }
    }

    Ibar.correctBoundaryConditions();

    return tIbar;
}


Foam::tmp<Foam::surfaceScalarField> Foam::localmultiscaleNeoHookeanElasticMisesPlastic::Ibar
(
    const surfaceSymmTensorField& devBEbar
)
{
    // From Simo & Hughes 1998:
    // but this incorrectly results in det(bEbar) =! 1
    //bEbar = (s/mu) + Ibar*I;

    // A method of calculating Ibar to enforce det(bEbar) == 1 is proposed
    // by solving a cubic equation.
    // Rubin and Attia, CALCULATION OF HYPERELASTIC RESPONSE OF FINITELY
    // DEFORMED ELASTIC-VISCOPLASTIC MATERIALS, INTERNATIONAL JOURNAL FOR
    // NUMERICAL METHODS IN ENGINEERING, VOL. 39,309-320(1996)
    // and
    // M. Hollenstein M. Jabareen M. B. Rubin, Modeling a smooth elastic-
    // inelastic transition with a strongly objective numerical integrator
    // needing no iteration, Comput Mech (2013) 52:649–667
    // DOI 10.1007/s00466-013-0838-7

    // Note: In Hollenstrain et al. (2013), they suggest that Eqn 49(a) in the
    // original Rubin and Attia paper should be used.

    // Method implemented below is translated from the SmoothMultiPhase fortran
    // subroutine of Rubin

    tmp<surfaceScalarField> tIbar
    (
        new surfaceScalarField
        (
            IOobject
            (
                "Ibar",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    surfaceScalarField& Ibar = tIbar();

    // Take reference to internal fields for efficiency
    scalarField& IbarI = Ibar.internalField();
    const symmTensorField devBEbarI = devBEbar.internalField();

    // Calculate internal field
    forAll(IbarI, cellI)
    {
        const scalar detdevBepr = det(devBEbarI[cellI]);
        const scalar dotprod = devBEbarI[cellI] && devBEbarI[cellI];
        const scalar fac1 = 2.0*dotprod/3.0;

        scalar alpha1 = 0.0;

        if (mag(fac1) < SMALL)
        {
            alpha1 = 3.0;
        }
        else
        {
            const scalar fac2 = (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

            if (fac2 >= 1.0)
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cosh(Foam::acosh(fac2)/3.0);
            }
            else
            {
                alpha1 = 3.0*Foam::sqrt(fac1)*Foam::cos(Foam::acos(fac2)/3.0);
            }
        }

        IbarI[cellI] = alpha1/3.0;
    }

    // Calculate boundary field
    forAll(Ibar.boundaryField(), patchI)
    {
        if
        (
            !Ibar.boundaryField()[patchI].coupled()
         && Ibar.boundaryField()[patchI].type() != "empty"
        )
        {
            // Take reference to patch fields for efficiency
            scalarField& IbarP = Ibar.boundaryField()[patchI];
            const symmTensorField& devBEbarP =
                devBEbar.boundaryField()[patchI];

            forAll(IbarP, faceI)
            {
                const scalar detdevBepr = det(devBEbarP[faceI]);
                const scalar dotprod =
                    devBEbarP[faceI] && devBEbarP[faceI];
                const scalar fac1 = 2.0*dotprod/3.0;

                scalar alpha1 = 0.0;

                if (mag(fac1) < SMALL)
                {
                    alpha1 = 3.0;
                }
                else
                {
                    const scalar fac2 =
                        (4.0*(1.0 - detdevBepr))/(pow(fac1, 1.5));

                    if (fac2 >= 1.0)
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cosh(Foam::acosh(fac2)/3.0);
                    }
                    else
                    {
                        alpha1 =
                            3.0*Foam::sqrt(fac1)
                           *Foam::cos(Foam::acos(fac2)/3.0);
                    }
                }

                IbarP[faceI] = alpha1/3.0;
            }
        }
    }

    Ibar.correctBoundaryConditions();

    return tIbar;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::localmultiscaleNeoHookeanElasticMisesPlastic::localmultiscaleNeoHookeanElasticMisesPlastic
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict,
    const nonLinearGeometry::nonLinearType& nonLinGeom
)
:
    pRveSolvers_(),
    applyBoundaryRve_
    (
           dict.lookupOrDefault<bool>("applyBoundaryRve", true)
    ),
    //RVEs_(1),
    mechanicalLaw(name, mesh, dict, nonLinGeom),
    rho_(dict.lookup("rho")),
    mu_("zero", dimPressure, 0.0),
    K_("zero", dimPressure, 0.0),
    relFPtr_(NULL),
    relFfPtr_(NULL),
    JPtr_(NULL),
    JfPtr_(NULL),
    FPtr_(NULL),
    FfPtr_(NULL),
    stressPlasticStrainSeries_(dict),
    sigmaHyd_
    (
        IOobject
        (
            "sigmaHyd",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimPressure, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    smoothPressure_(dict.lookupOrDefault<Switch>("smoothPressure", true)),
    sigmaY_
    (
        IOobject
        (
            "sigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "initialYieldStress", dimPressure, stressPlasticStrainSeries_(0.0)
        )
    ),
    sigmaYf_
    (
        IOobject
        (
            "sigmaYf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
           "initialYieldStress", dimPressure, stressPlasticStrainSeries_(0.0)
        )
    ),
    DSigmaY_
    (
        IOobject
        (
            "DSigmaY",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    DSigmaYf_
    (
        IOobject
        (
            "DSigmaYf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimPressure, 0.0)
    ),
    epsilonP_
    (
        IOobject
        (
            "epsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    epsilonPf_
    (
        IOobject
        (
            "epsilonPf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonP_
    (
        IOobject
        (
            "DEpsilonP",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    DEpsilonPf_
    (
        IOobject
        (
            "DEpsilonPf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    bEbarTrial_
    (
        IOobject
        (
            "bEbarTrial",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbarTrialf_
    (
        IOobject
        (
            "bEbarTrialf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbar_
    (
        IOobject
        (
            "bEbar",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    bEbarf_
    (
        IOobject
        (
            "bEbarf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("I", dimless, I)
    ),
    DEpsilonPEq_
    (
        IOobject
        (
            "DEpsilonPEq",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    DEpsilonPEqf_
    (
        IOobject
        (
            "DEpsilonPEqf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    DLambda_
    (
        IOobject
        (
            "DLambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    DLambdaf_
    (
        IOobject
        (
            "DLambdaf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    epsilonPEq_
    (
        IOobject
        (
            "epsilonPEq",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    epsilonPEqf_
    (
        IOobject
        (
            "epsilonPEqf",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0.0)
    ),
    activeYield_
    (
        IOobject
        (
            "activeYield",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless, 0)
    ),
    plasticN_
    (
        IOobject
        (
            "plasticN",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    plasticNf_
    (
        IOobject
        (
            "plasticNf",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    ),
    nonLinearPlasticity_(stressPlasticStrainSeries_.size() > 2),
    updateBEbarConsistent_
    (
        dict.lookupOrDefault<Switch>
        (
            "updateBEbarConsistent",
            Switch(true)
        )
    ),
    Hp_(0.0),
    maxDeltaErr_
    (
        mesh.time().controlDict().lookupOrDefault<scalar>("maxDeltaErr", 0.01)
    )
{
    // Force storage of old time for adjustable time-step calculations
    plasticN_.oldTime();
    bEbar_.oldTime();
    
    // Determine number of rve (cells + non-empty boundary faces)
    label nRve = mesh.nCells();

    Pout << "applyBoundaryRve: " <<  applyBoundaryRve_ << endl;

    if (applyBoundaryRve_)
    {
    forAll(mesh.boundary(), patchI)
    {
        //if (!mesh.boundary()[patchI].coupled())
        if (mesh.boundaryMesh()[patchI].name() == "floor")
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

    Info<< "    smoothPressure: " << smoothPressure_ << nl
        << "    updateBEbarConsistent: " << updateBEbarConsistent_ << endl;

    // Read elastic parameters
    // The user can specify E and nu or mu and K
    if (dict.found("E") && dict.found("nu"))
    {
        // Read the Young's modulus
        const dimensionedScalar E = dimensionedScalar(dict.lookup("E"));

        // Read the Poisson's ratio
        const dimensionedScalar nu = dimensionedScalar(dict.lookup("nu"));

        // Set the shear modulus
        mu_ = E/(2.0*(1.0 + nu));

        // Set the bulk modulus
        if (planeStress())
        {
            K_ = (nu*E/((1.0 + nu)*(1.0 - nu))) + (2.0/3.0)*mu_;
        }
        else
        {
            K_ = (nu*E/((1.0 + nu)*(1.0 - 2.0*nu))) + (2.0/3.0)*mu_;
        }
    }
    else if (dict.found("mu") && dict.found("K"))
    {
        mu_ = dimensionedScalar(dict.lookup("mu"));
        K_ = dimensionedScalar(dict.lookup("K"));
    }
    else
    {
        FatalErrorIn
        (
            "multiscalNeoHookeanElasticMisesPlastic::multiscalNeoHookeanElasticMisesPlastic::()"
        )   << "Either E and nu or mu and K elastic parameters should be "
            << "specified" << abort(FatalError);
    }

    // Check if plasticity is a nonlinear function of plastic strain
    if (nonLinearPlasticity_)
    {
        Info<< "    Plasticity is nonlinear" << endl;
    }
    else
    {
        if (stressPlasticStrainSeries_.size() == 1)
        {
            Info<< "    Perfect Plasticity" << endl;
        }
        else
        {
            Info<< "    Plasticity is linear" << endl;

            // Define linear plastic modulus
            Hp_ =
                (
                    stressPlasticStrainSeries_[1].second()
                  - stressPlasticStrainSeries_[0].second()
                )
               /(
                    stressPlasticStrainSeries_[1].first()
                  - stressPlasticStrainSeries_[0].first()
                );
        }
    }

    if (updateBEbarConsistent_)
    {
        Info<< "updateBEbarConsistent is active" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::localmultiscaleNeoHookeanElasticMisesPlastic::~localmultiscaleNeoHookeanElasticMisesPlastic()
{
    deleteDemandDrivenData(relFPtr_);
    deleteDemandDrivenData(relFfPtr_);
    deleteDemandDrivenData(JPtr_);
    deleteDemandDrivenData(JfPtr_);
    deleteDemandDrivenData(FPtr_);
    deleteDemandDrivenData(FfPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::localmultiscaleNeoHookeanElasticMisesPlastic::rho() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "lawRho",
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


Foam::tmp<Foam::volScalarField>
Foam::localmultiscaleNeoHookeanElasticMisesPlastic::impK() const
{
    // Calculate scaling factor to ensure optimal convergence
    // This is similar to the tangent matrix in FE procedures

    // Calculate deviatoric trial stress
    const volSymmTensorField sTrial = mu_*dev(bEbarTrial_);

    const volScalarField Ibar = tr(bEbarTrial_)/3.0;
    const volScalarField muBar = Ibar*mu_;

    // Magnitude of the deviatoric trial stress
    const volScalarField magSTrial =
        max(mag(sTrial), dimensionedScalar("SMALL", dimPressure, SMALL));

    // Calculate scaling factor
    const volScalarField scaleFactor = 1.0 - (2.0*muBar*DLambda_/magSTrial);

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
            //mesh(),
            //(4.0/3.0)*mu_ + K_, // == 2*mu + lambda
            //zeroGradientFvPatchScalarField::typeName
            scaleFactor*(4.0/3.0)*mu_ + K_
        )
    );
}


void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::correct(volSymmTensorField& sigma)
{
    // Check if the mathematical model is in total or updated Lagrangian form
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        if (!incremental())
        {
            FatalErrorIn(type() + "::correct(volSymmTensorField& sigma)")
                << "Not implemented for non-incremental updated Lagrangian"
                << abort(FatalError);
        }

        // Lookup gradient of displacement increment
        const volTensorField& gradDD =
            mesh().lookupObject<volTensorField>("grad(DD)");

        // Update the relative deformation gradient
        relF() = I + gradDD.T();

        // Update the total deformation gradient
        F() = relF() & F().oldTime();
    }
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
            // Lookup gradient of displacement
            const volTensorField& gradD =
                mesh().lookupObject<volTensorField>("grad(D)");

            // Update the total deformation gradient
            F() = I + gradD.T();

            // Update the relative deformation gradient
            relF() = F() & inv(F().oldTime());
            
            forAll(F().internalField(), cellI)
            {
                pRveSolvers_[cellI].avgRelF() = relF().internalField()[cellI];
                //Info << "reF() is : " << relF().internalField()[cellI] << endl;
            }
        
            if (applyBoundaryRve_)
            {
            //label solI = mesh().nCells();
            label nRve = mesh().nCells();
            //Pout << "Number of RVEs 1 : " << nRve << endl;
            forAll(relF().boundaryField(), patchI)
              {
                forAll(relF().boundaryField()[patchI],faceI)
                {
                        pRveSolvers_[nRve].avgRelF() =
                          relF().boundaryField()[patchI][faceI];
                nRve++;
                }
              }
            }
    }
    else
    {
        FatalErrorIn
        (
            type() + "::correct(volSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }

    // Update the Jacobian of the total deformation gradient
    J() = det(F());

    // Calculate the relative Jacobian
    const volScalarField relJ = J()/J().oldTime();

    // Calculate the relative deformation gradient with the volumetric term
    // removed
    const volTensorField relFbar = pow(relJ, -1.0/3.0)*relF();

    // Update bE trial
    bEbarTrial_ = transform(relFbar, bEbar_.oldTime());

    // Calculate trial deviatoric stress
    const volSymmTensorField sTrial = mu_*dev(bEbarTrial_);

    const volScalarField Ibar = tr(bEbarTrial_)/3.0;
    const volScalarField muBar = Ibar*mu_;

    // Check for plastic loading
    // and calculate increment of plastic equivalent strain
    // i.e. the plastic multiplier

    // Store previous iteration for under-relaxation and calculation of plastic
    // residual in the solver
    DEpsilonP_.storePrevIter();

    // Normalise residual in Newton method with respect to mag(bE)
    const scalar maxMagBE = max(gMax(mag(bEbarTrial_.internalField())), SMALL);

    // Trial yield function
    // sigmaY is the Cauchy yield stress so we scale it by J
    const volScalarField fTrial = mag(sTrial) - sqrtTwoOverThree_*J()*sigmaY_;

    // Magnitude of hardening slope
    const scalar magHp = mag(Hp_);

    // Take references to the internal fields for efficiency
    const scalarField& fTrialI = fTrial.internalField();
    const symmTensorField& sTrialI = sTrial.internalField();
    symmTensorField& plasticNI = plasticN_.internalField();
    scalarField& DSigmaYI = DSigmaY_.internalField();
    scalarField& DLambdaI = DLambda_.internalField();
    const scalarField& muBarI = muBar.internalField();
    const scalarField& JI = J().internalField();
    const scalarField& sigmaYI = sigmaY_.internalField();
    const scalarField& epsilonPEqOldI = epsilonPEq_.oldTime().internalField();

    // Calculate DLambda_ and plasticN_
    forAll(fTrialI, cellI)
    {
        // Calculate return direction plasticN
        const scalar magS = mag(sTrialI[cellI]);
        if (magS > SMALL)
        {
            plasticNI[cellI] = sTrialI[cellI]/magS;
        }

        // Calculate DLambda/DEpsilonPEq
        if (fTrialI[cellI] < SMALL)
        {
            // elastic
            DSigmaYI[cellI] = 0.0;
            DLambdaI[cellI] = 0.0;
        }
        else
        {
            if (nonLinearPlasticity_)
            {
                // Total equivalent plastic strain where t is start of time-step
                scalar curSigmaY = 0.0; // updated in loop below

                // Calculates DEpsilonPEq using Newtons's method
                newtonLoop
                (
                    DLambdaI[cellI],
                    curSigmaY,
                    epsilonPEqOldI[cellI],
                    magS,
                    muBarI[cellI],
                    JI[cellI],
                    maxMagBE
                );

                // Update increment of yield stress
                DSigmaYI[cellI] = curSigmaY - sigmaYI[cellI];
            }
            else
            {
                // Plastic modulus is linear
                DLambdaI[cellI] = fTrialI[cellI]/(2*muBarI[cellI]);

                if (magHp > SMALL)
                {
                    DLambdaI[cellI] /= 1.0 + Hp_/(3*muBarI[cellI]);

                    // Update increment of yield stress
                    DSigmaYI[cellI] = DLambdaI[cellI]*Hp_;
                }
            }
        }
    }

    forAll(fTrial.boundaryField(), patchI)
    {
        // Take references to the boundary patch fields for efficiency
        const scalarField& fTrialP = fTrial.boundaryField()[patchI];
        const symmTensorField& sTrialP = sTrial.boundaryField()[patchI];
        symmTensorField& plasticNP = plasticN_.boundaryField()[patchI];
        scalarField& DSigmaYP = DSigmaY_.boundaryField()[patchI];
        scalarField& DLambdaP = DLambda_.boundaryField()[patchI];
        const scalarField& muBarP = muBar.boundaryField()[patchI];
        const scalarField& JP = J().boundaryField()[patchI];
        const scalarField& sigmaYP = sigmaY_.boundaryField()[patchI];
        const scalarField& epsilonPEqOldP =
            epsilonPEq_.oldTime().boundaryField()[patchI];

        forAll(fTrialP, faceI)
        {
            // Calculate direction plasticN
            const scalar magS = mag(sTrialP[faceI]);
            if (magS > SMALL)
            {
                plasticNP[faceI] = sTrialP[faceI]/magS;
            }

            // Calculate DEpsilonPEq
            if (fTrialP[faceI] < SMALL)
            {
                // elasticity
                DSigmaYP[faceI] = 0.0;
                DLambdaP[faceI] = 0.0;
            }
            else
            {
                // yielding
                if (nonLinearPlasticity_)
                {
                    scalar curSigmaY = 0.0; // updated in loop below

                    // Calculate DEpsilonPEq and curSigmaY
                    newtonLoop
                    (
                        DLambdaP[faceI],
                        curSigmaY,
                        epsilonPEqOldP[faceI],
                        magS,
                        muBarP[faceI],
                        JP[faceI],
                        maxMagBE
                    );

                    // Update increment of yield stress
                    DSigmaYP[faceI] = curSigmaY - sigmaYP[faceI];
                }
                else
                {
                    // Plastic modulus is linear
                    DLambdaP[faceI] = fTrialP[faceI]/(2.0*muBarP[faceI]);

                    if (magHp > SMALL)
                    {
                        DLambdaP[faceI] /= 1.0 + Hp_/(3.0*muBarP[faceI]);

                        // Update increment of yield stress
                        DSigmaYP[faceI] = DLambdaP[faceI]*Hp_;
                    }
                }
            }
        }
    }

    // Update DEpsilonP and DEpsilonPEq
    DEpsilonPEq_ = sqrtTwoOverThree_*DLambda_;
    DEpsilonP_ = Ibar*DLambda_*plasticN_;
    DEpsilonP_.relax();

    // Calculate deviatoric stress
    const volSymmTensorField s = sTrial - 2*mu_*DEpsilonP_;

    // Update bEbar
    if (updateBEbarConsistent_)
    {
        const volSymmTensorField devBEbar = (s/mu_);
        bEbar_ = devBEbar + this->Ibar(devBEbar)*I;
    }
    else
    {
        bEbar_ = (s/mu_) + Ibar*I;
    }

    // Update hydrostatic stress (negative of pressure)
    if (smoothPressure_)
    {
        // Calculate the hydrostatic pressure by solving a Laplace equation;
        // this ensures smoothness of the field and quells oscillations

        // Lookup the momentum equation inverse diagonal field
        const volScalarField AD = mesh().lookupObject<volScalarField>("DEqnA");

        // Pressure diffusivity field
        // Note: (4.0/3.0)*mu + K == 2*mu + lambda
        const surfaceScalarField rDAf
        (
            "rDAf",
            fvc::interpolate
            (
                ((4.0/3.0)*mu_ + K_)/AD, "interpolate(grad(sigmaHyd))"
            )
        );

        const dimensionedScalar one("one", dimless, 1.0);
        //const dimensionedScalar fac(dict().lookup("smoothFactor"));

        // Construct the pressure equation
        fvScalarMatrix sigmaHydEqn
        (
            fvm::Sp(one, sigmaHyd_)
          - fvm::laplacian(rDAf, sigmaHyd_, "laplacian(DDD,DD)")
         ==
            0.5*K_*(pow(J(), 2.0) - 1.0)
          - fvc::div(rDAf*fvc::interpolate(fvc::grad(sigmaHyd_)) & mesh().Sf())
        );

        // Store the pressue field to allow under-relaxation
        sigmaHyd_.storePrevIter();

        // Under-relax the linear system
        sigmaHydEqn.relax(0.7);

        // Solve the pressure equation
        sigmaHydEqn.solve();

        // Under-relax the pressure field
        sigmaHyd_.relax(0.2);
    }
    else
    {
        // Calculate the hydrostatic pressure directly from the displacement
        // field
        sigmaHyd_ = 0.5*K_*(pow(J(), 2.0) - 1.0);
    }

    // Update the Cauchy stress
    sigma = (1.0/J())*(sigmaHyd_*I + s);

    //evelve the RVE with plasticity cell

    //label iRVEs = 0;

    forAll(pRveSolvers_, nRve)
    {
    	if(activeYield_.internalField()[nRve] == 1.0)
    	{
    		Info << "Solving pRve: " << nRve << endl;
                pRveSolvers_[nRve].evolve();
    	}
    }
    //Info << "Number of RVEs applied in this time: "<< iRVEs << endl;
    
    //forAll(pRveSolvers_, nRve)
    //{
    //    Info << "Solving pRve: " << nRve << endl;
    //    pRveSolvers_[nRve].evolve();
    //}
    
    // Get average stess from RVE model
    

    forAll(sigma.internalField(), cellI)
    {
    	if(activeYield_.internalField()[cellI] == 1.0)
    	{
    		Info << "Get average stess from RVE: " << cellI << endl;
                sigma.internalField()[cellI] = pRveSolvers_[cellI].avgSigma().value();
    	}
    }

    //RVEs_.clear();
    

    //label iRVEs = 0;
    //forAll(C.internalField(), cellI)
    //{
    //	if(activeYield_.internalField()[cellI] == 1.0)
    //	{
    //		RVE_[iRVEs][0] = cellI;
	//	RVE_[iRVEs][1] = C.internalField()[cellI][0];
     //           RVE_[iRVEs][2] = C.internalField()[cellI][1];
     //           RVE_[iRVEs][3] = C.internalField()[cellI][2];
     //           iRVEs++;
    //	}
    //}

    //forAll(sigma.internalField(), cellI)
    //{
    //            sigma.internalField()[cellI] =
    //            pRveSolvers_[cellI].avgSigma().value();
    //}
    
    if (applyBoundaryRve_)
    {
    //solI = mesh().nCells();
    label nRve = mesh().nCells();
    forAll(sigma.boundaryField(), patchI)
      {
            forAll(sigma.boundaryField()[patchI], faceI)
            {
                sigma.boundaryField()[patchI][faceI] =
                    pRveSolvers_[nRve].avgSigma().value();
                nRve++;
                //solI++;
            }
      }
    }
    //sigma.correctBoundaryConditions();
}


void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::correct(surfaceSymmTensorField& sigma)
{
    // Check if the mathematical model is in total or updated Lagrangian form
    if (nonLinGeom() == nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
        if (!incremental())
        {
            FatalErrorIn(type() + "::correct(surfaceSymmTensorField& sigma)")
                << "Not implemented for non-incremental updated Lagrangian"
                << abort(FatalError);
        }

        // Lookup gradient of displacement increment
        const surfaceTensorField& gradDD =
            mesh().lookupObject<surfaceTensorField>("grad(DD)f");

        // Update the relative deformation gradient
        relFf() = I + gradDD.T();

        // Update the total deformation gradient
        Ff() = relFf() & Ff().oldTime();
    }
    else if (nonLinGeom() == nonLinearGeometry::TOTAL_LAGRANGIAN)
    {
            // Lookup gradient of displacement
            const surfaceTensorField& gradD =
                mesh().lookupObject<surfaceTensorField>("grad(D)f");

            // Update the total deformation gradient
            Ff() = I + gradD.T();

            // Update the relative deformation gradient
            relFf() = Ff() & inv(Ff().oldTime());
    }
    else
    {
        FatalErrorIn
        (
            type() + "::correct(surfaceSymmTensorField& sigma)"
        )   << "Unknown nonLinGeom type: " << nonLinGeom() << abort(FatalError);
    }
}


Foam::scalar Foam::localmultiscaleNeoHookeanElasticMisesPlastic::residual()
{
    // Calculate residual based on change in plastic strain increment
    if
    (
        mesh().db().lookupObject<fvMesh>
        (
            baseMeshRegionName()
        ).foundObject<surfaceTensorField>("Ff")
    )
    {
        return
            gMax
            (
                mag
                (
                    DEpsilonPf_.internalField()
                  - DEpsilonPf_.prevIter().internalField()
                )
            )/gMax(SMALL + mag(DEpsilonPf_.prevIter().internalField()));
    }
    else
    {
        return
            gMax
            (
                mag
                (
                    DEpsilonP_.internalField()
                  - DEpsilonP_.prevIter().internalField()
                )
            )/gMax(SMALL + mag(DEpsilonP_.prevIter().internalField()));
    }
}


void Foam::localmultiscaleNeoHookeanElasticMisesPlastic::updateTotalFields()
{
    Info<< nl << "Updating total accumulated fields" << endl;
    sigmaY_ += DSigmaY_;
    sigmaYf_ += DSigmaYf_;

    Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;
    epsilonPEq_ += DEpsilonPEq_;
    epsilonPEqf_ += DEpsilonPEqf_;
    epsilonP_ += DEpsilonP_;
    epsilonPf_ += DEpsilonPf_;

    // Count cells actively yielding
    int numCellsYielding = 0;

    forAll(activeYield_.internalField(), celli)
    {
        if (DEpsilonPEq_.internalField()[celli] > SMALL)
        {
            activeYield_.internalField()[celli] = 1.0;
            numCellsYielding++;
        }
        else
        {
            activeYield_.internalField()[celli] = 0.0;
        }
    }

    const volVectorField& C = mesh().C();
    forAll(C.internalField(), cellI)
    {
    	if(activeYield_.internalField()[cellI] == 1.0)
    	{
    		Info << "The center of RVE-: " << cellI <<"is"<< C.internalField()[cellI]<< endl;
    	}
    }

    reduce(numCellsYielding, sumOp<int>());

    forAll(activeYield_.boundaryField(), patchi)
    {
        if (!activeYield_.boundaryField()[patchi].coupled())
        {
            forAll(activeYield_.boundaryField()[patchi], facei)
            {
                if (DEpsilonPEq_.boundaryField()[patchi][facei] > SMALL)
                {
                    activeYield_.boundaryField()[patchi][facei] = 1.0;
                }
                else
                {
                    activeYield_.boundaryField()[patchi][facei] = 0.0;
                }
            }
        }
    }

    activeYield_.correctBoundaryConditions();

    Info<< "    " << numCellsYielding << " cells are actively yielding"
        << nl << endl;
    
    forAll(pRveSolvers_, nRve)
    {
        pRveSolvers_[nRve].updateTotalFields();
    }
}


Foam::scalar Foam::localmultiscaleNeoHookeanElasticMisesPlastic::newDeltaT()
{
    // In the calculation of the plastic strain increment, the return direction
    // is kept constant for the time-step; we can approximate the error based on
    // the difference in the return direction from the start to the end of the
    // time-step, where the return direction is given normalised deviatoric
    // strain. The error approximation is obtained using the difference between
    // the trapezoidal rule and the EUler backward method, as described in:

    // Nam-Sua Lee, Klaus-Jurgen Bathe, Error indicators and adaptive remeshing
    // in large deformation finite element analysis, Finite Elements in
    // Analysis and Design 16 (1994) 99-139.

    // Update the total deformatio gradient
    if (mesh().foundObject<surfaceTensorField>("grad(DD)f"))
    {
        F() = fvc::average(relFf()) & F().oldTime();
    }
    else
    {
        F() = relF() & F().oldTime();
    }

    // Calculate the total true (Hencky) strain
    const volSymmTensorField epsilon = 0.5*log(symm(F().T() & F()));

    // Calculate equivalent strain, for normalisation of the error
    const volScalarField epsilonEq = sqrt((2.0/3.0)*magSqr(dev(epsilon)));

    // Take reference to internal fields
    const symmTensorField& DEpsilonPI = DEpsilonP_.internalField();
    const symmTensorField& plasticNI = plasticN_.internalField();
    const symmTensorField& plasticNIold = plasticN_.oldTime().internalField();
    const scalarField& epsilonEqI = epsilonEq.internalField();

    // Calculate error field
    const symmTensorField DEpsilonPErrorI =
        Foam::sqrt(3.0/8.0)*DEpsilonPI*mag(plasticNI - plasticNIold)
       /(epsilonEqI + SMALL);

    // Max error
    const scalar maxMagDEpsilonPErr = gMax(mag(DEpsilonPErrorI));

    if (maxMagDEpsilonPErr > SMALL)
    {
        Info<< "    " << name() << ": max time integration error = "
            << maxMagDEpsilonPErr
            << endl;

        if (maxMagDEpsilonPErr > 50*maxDeltaErr_)
        {
            WarningIn
            (
                "Foam::scalar Foam::neoHookeanElasticMisesPlastic::newDeltaT()"
                " const"
            )   << "The error in the plastic strain is lover 50 times larger "
                << "than the desired value!\n    Consider starting the "
                << "simulation with a smaller initial time-step" << endl;
        }

        // Calculate the time-step scaling factor, where maxDeltaErr_ is the
        // maximum allowed error
        const scalar scaleFac = maxDeltaErr_/maxMagDEpsilonPErr;

        // Return the new time-step size
        return scaleFac*mesh().time().deltaTValue();
    }

    return mesh().time().endTime().value();
}


// ************************************************************************* //
