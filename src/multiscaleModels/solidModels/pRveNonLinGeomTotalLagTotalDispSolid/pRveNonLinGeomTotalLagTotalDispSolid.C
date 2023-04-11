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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "pRveNonLinGeomTotalLagTotalDispSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"

#include "skewCorrectionVectors.H"
#include "gaussGrad.H"
#include "logVolFields.H"
#include "polyPatchID.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pRveNonLinGeomTotalLagTotalDispSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, pRveNonLinGeomTotalLagTotalDispSolid, dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //


void pRveNonLinGeomTotalLagTotalDispSolid::predict()
{
    Info << "Predicting solid model" << endl;

    // Predict D using the velocity field
    // Note: the case may be steady-state but U can still be calculated using a
    // transient method
    D() = D().oldTime() + U()*runTime().deltaT();

    // Update gradient of displacement
    mechanical().grad(D(), gradD());

    // Total deformation gradient
    F_ = I + gradD().T();

    // Inverse of the deformation gradient
    Finv_ = inv(F_);

    // // Jacobian of the deformation gradient
    J_ = det(F_);

    // Calculate the stress using run-time selectable mechanical law
    mechanical().correct(sigma());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pRveNonLinGeomTotalLagTotalDispSolid::pRveNonLinGeomTotalLagTotalDispSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    avgFseries_(),
    avgF_
    (
        IOobject
        (
            "avgF",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor
        (
            "I",
            dimless,
            // I
            tensor(1, 0, 0, -2.75331e-05, 1.00003, -2.23866e-22, 0, 0, 1)
        )
    ),
    avgSigma_("0", dimPressure, symmTensor::zero),
    totD_
    (
        IOobject
        (
            "totD",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimLength, vector::zero)
    ),
    pointTotD_
    (
        IOobject
        (
            "pointTotD",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimLength, vector::zero)
    ),
    avgD_
    (
        IOobject
        (
            "avgD",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimLength, vector::zero)
    ),
    F_
    (
        IOobject
        (
            "F",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedTensor("I", dimless, I)
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(F_)
    ),
    J_
    (
        IOobject
        (
            "J",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        det(F_)
    ),
    impK_(mechanical().impK()),
    impKf_(mechanical().impKf()),
    rImpK_(1.0/impK_),
    df_
    (
        IOobject
        (
            "df",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh(),
        dimensionedVector("0", dimLength, vector::zero)
    ),
    predict_(solidProperties().lookupOrDefault<Switch>("predict", false))
{
    avgD_.oldTime();

    // Check if average deformation gradient is time-varying
    if (solidProperties().found("avgDeformationGradientSeries"))
    {
        Info<< "    average deformation gradient is time-varying" << endl;
        avgFseries_ =
            interpolationTable<tensor>
            (
                solidProperties().subDict("avgDeformationGradientSeries")
            );

        avgF_ = avgFseries_(runTime.value());
    }

    // Calculate df
    const labelList& own = mesh().owner();
    const labelList& nei = mesh().neighbour();
    const vectorField& C = mesh().C().internalField();
    forAll(df_.internalField(), faceI)
    {
        df_.internalField()[faceI] = (C[nei[faceI]] - C[own[faceI]]);
    }
    forAll(df_.boundaryField(), patchI)
    {
        forAll(df_.boundaryField()[patchI], faceI)
        {
            df_.boundaryField()[patchI] =
                mesh().boundary()[patchI].nf()
               /mesh().deltaCoeffs().boundaryField()[patchI];
        }
    }      
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool pRveNonLinGeomTotalLagTotalDispSolid::evolve()
{
    Info<< "Evolving solid solver: " << this->type() << endl;

    if (predict_)
    {
        predict();
    }

    int iCorr = 0;
    lduSolverPerformance solverPerfD;
    blockLduMatrix::debug = 0;

    Info<< "Solving the total Lagrangian form of the momentum equation for D"
        << endl;

    if (avgFseries_.size())
    {
        avgF_ = avgFseries_(runTime().value());
    }

    Info << "Avg deformation gradient: "
         << average(avgF_.internalField()) << endl;
        
    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        D().storePrevIter();

        avgD_ = ((avgF_ - I) & mesh().C());
        avgD_.correctBoundaryConditions();
        
        // Momentum equation total displacement total Lagrangian form
        fvVectorMatrix DEqn
        (
            rho()*fvm::d2dt2(D())
          + rho()*fvc::d2dt2(avgD_)
         == fvm::laplacian(impKf_, D(), "laplacian(DD,D)")
          - fvc::laplacian(impKf_, D(), "laplacian(DD,D)")
          + fvc::div(J_*Finv_ & sigma(), "div(sigma)")
          + mechanical().RhieChowCorrection(D(), gradD())
        );

        // Under-relax the linear system
        DEqn.relax();

        // Solve the linear system
        solverPerfD = DEqn.solve();

        // Fixed or adaptive field under-relaxation
        // relaxField(D(), iCorr);
        D().relax();

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Update gradient of displacement
        mechanical().grad(D(), gradD());

        // Update gradient of displacement increment
        gradDD() = gradD() - gradD().oldTime();

        // Total deformation gradient
        F_ = avgF_ + gradD().T();

        // Inverse of the deformation gradient
        Finv_ = inv(F_);

        // Jacobian of the deformation gradient
        J_ = det(F_);

        // Update the momentum equation inverse diagonal field
        // This may be used by the mechanical law when calculating the
        // hydrostatic pressure
        const volScalarField DEqnA("DEqnA", DEqn.A());

        // Add average gradient to gradient of displacement perturbation
        gradD() += (avgF_.T() - I);

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());

        // Return to gradient of displacement perturbation
        gradD() -= (avgF_.T() - I);
    }
    while
    (
       !converged
        (
            iCorr,
            solverPerfD.initialResidual(),
            solverPerfD.nIterations(),
            D()
        ) && ++iCorr < nCorr()
    );
    
    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), pointD());

    pointTotD_.internalField() = ((avgF_[0] - I) & mesh().points()) + pointD();

    // Increment of point displacement
    pointDD() = pointD() - pointD().oldTime();

    // Velocity
    U() = fvc::ddt(D());

    // Update average Cauchy stress
    {
        scalar totalVolume = 0;

        forAll(mesh().boundary(), patchI)
        {
            if (isA<cyclicFvPatch>(mesh().boundary()[patchI]))
            {
                totalVolume +=
                    gSum
                    (
                        mesh().Cf().boundaryField()[patchI]
                      & mesh().Sf().boundaryField()[patchI]
                    )/2;
            }
        }
        
        // polyPatchID voidPatch(word("void"), mesh().boundaryMesh());

        // if (voidPatch.active())
        // {
        //     Vv =
        //       - gSum
        //         (
        //             mesh().Cf().boundaryField()[voidPatch.index()]
        //           & mesh().Sf().boundaryField()[voidPatch.index()]
        //         )/2;
        // }

        scalar Vs = gSum(mesh().V().field());

        Info << "Vs = " << Vs << ", totalVolume = " << totalVolume << endl;

        avgSigma_.value() =
            gSum(sigma().internalField()*mesh().V().field())
           /totalVolume;

        Info << "Avg Cauchy stress: "
             << avgSigma_.value() << endl;
    }
    
    return true;
}


tmp<vectorField> pRveNonLinGeomTotalLagTotalDispSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    // Patch index
    const label patchID = patch.index();

    // Patch implicit stiffness field
    const scalarField& impK = impK_.boundaryField()[patchID];

    // Patch reciprocal implicit stiffness field
    const scalarField& rImpK = rImpK_.boundaryField()[patchID];

    // Patch gradient
    const tensorField& pGradD = gradD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch total deformation gradient inverse
    const tensorField& Finv = Finv_.boundaryField()[patchID];

    // Patch Jacobian
    const scalarField& J = J_.boundaryField()[patchID];
    
    // Patch unit normals (initial configuration)
    const vectorField n = patch.nf();

    // Patch unit normals (deformed configuration)
    vectorField nCurrent = Finv.T() & n;
    nCurrent /= mag(nCurrent);

    // Return patch snGrad
    return tmp<vectorField>
    (
        new vectorField
        (
            (
                (traction - nCurrent*pressure)
              - (n & (Finv & (J*pSigma)))
              // - (nCurrent & pSigma)
              + impK*(n & pGradD)
            )*rImpK
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
