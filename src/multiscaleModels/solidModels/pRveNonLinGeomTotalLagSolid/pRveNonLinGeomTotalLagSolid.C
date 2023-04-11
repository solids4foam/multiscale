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

#include "pRveNonLinGeomTotalLagSolid.H"
#include "volFields.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "solidTractionFvPatchVectorField.H"
#include "fvcGradf.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pRveNonLinGeomTotalLagSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, pRveNonLinGeomTotalLagSolid, dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pRveNonLinGeomTotalLagSolid::pRveNonLinGeomTotalLagSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    avgRelF_
    (
        "I",
        dimless,
        I
        // tensor(1, 0, 0, -2.75331e-05, 1.00003, -2.23866e-22, 0, 0, 1)
    ),
    avgF_
    (
        "I",
        dimless,
        I
        // tensor(1, 0, 0, -2.75331e-05, 1.00003, -2.23866e-22, 0, 0, 1)
    ),
    avgP_("0", dimPressure, tensor::zero),
    avgSigma_("0", dimPressure, symmTensor::zero),
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
            IOobject::NO_WRITE
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
    timeIndex_(-1)
{
    gradD().oldTime();

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


bool pRveNonLinGeomTotalLagSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    int iCorr = 0;
    lduSolverPerformance solverPerfDD;
    blockLduMatrix::debug = 0;

    Info<< "Solving the total Lagrangian form of the momentum equation for DD"
        << endl;

    // Reset enforceLinear switch
    enforceLinear() = false;

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        DD().storePrevIter();

        // Momentum equation incremental displacement total Lagrangian form
        fvVectorMatrix DDEqn
        (
            fvm::d2dt2(rho(), DD())
          + fvc::d2dt2(rho().oldTime(), D().oldTime())
         == fvm::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          - fvc::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          + fvc::div(J_*Finv_ & sigma(), "div(sigma)")
          + rho()*g()
          + mechanical().RhieChowCorrection(DD(), gradDD())
        );

        // Under-relax the linear system
        DDEqn.relax();

        // Solve the linear system
        solverPerfDD = DDEqn.solve();

        // Under-relax the DD field using fixed or adaptive under-relaxation
        relaxField(DD(), iCorr);

        // Update the total displacement
        D() = D().oldTime() + DD();

        // Update gradient of displacement increment
        mechanical().grad(DD(), gradDD());

        // Add average part
        gradDD() += (avgRelF_.value().T() - I);

        // Update the gradient of total displacement
        gradD() = gradD().oldTime() + gradDD();

        // // Total deformation gradient
        // F_ = relF & F_.oldTime();

        // Total deformation gradient
        F_ = I + gradD().T();

        // Inverse of the deformation gradient
        Finv_ = inv(F_);

        // Jacobian of the deformation gradient
        J_ = det(F_);

        // Update the momentum equation inverse diagonal field
        // This may be used by the mechanical law when calculating the
        // hydrostatic pressure
        const volScalarField DEqnA("DEqnA", DDEqn.A());

        // Calculate the stress using run-time selectable mechanical law
        mechanical().correct(sigma());
    }
    
    while
    (
       !converged
        (
            iCorr,
            solverPerfDD.initialResidual(),
            solverPerfDD.nIterations(),
            DD()
        ) && ++iCorr < nCorr()
    );

    // Interpolate cell displacements to vertices
    mechanical().interpolate(D(), pointD());

    // Increment of displacement
    D() = D().oldTime() + DD();

    // Increment of point displacement
    pointD() = pointD().oldTime() + pointDD();

    // Velocity
    U() = fvc::ddt(D());

    // Update average Cauchy stress
    updateAvgSigma(true);
    // {
    //     boundBox box(mesh().points());
        
    //     scalar totalVolume = box.span().x()*box.span().y()*box.span().z();

    //     // forAll(mesh().boundary(), patchI)
    //     // {
    //     //     if (isA<cyclicFvPatch>(mesh().boundary()[patchI]))
    //     //     {
    //     //         totalVolume +=
    //     //             gSum
    //     //             (
    //     //                 mesh().Cf().boundaryField()[patchI]
    //     //               & mesh().Sf().boundaryField()[patchI]
    //     //             )/2;
    //     //     }
    //     // }
        
    //     scalar Vs = gSum(mesh().V().field());

    //     Info << "Vs = " << Vs << ", totalVolume = " << totalVolume << endl;

    //     avgSigma_.value() =
    //         gSum(sigma().internalField()*mesh().V().field())
    //        /(totalVolume + SMALL);

    //     Info << "Avg Cauchy stress: "
    //          << avgSigma_.value() << endl;
    // }
    
    return true;
}
    
tmp<vectorField> pRveNonLinGeomTotalLagSolid::tractionBoundarySnGrad
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
    tensorField pGradDD = gradDD().boundaryField()[patchID];
    pGradDD -= (avgRelF_.value().T() - I);

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch unit normals (initial configuration)
    const vectorField n = patch.nf();

    // if (enforceLinear())
    // {
    //     // Return patch snGrad
    //     return tmp<vectorField>
    //     (
    //         new vectorField
    //         (
    //             (
    //                 (traction - n*pressure)
    //               - (n & pSigma)
    //               + impK*(n & pGradDD)
    //             )*rImpK
    //         )
    //     );
    // }
    // else
    {
        // Patch total deformation gradient inverse
        const tensorField& Finv = Finv_.boundaryField()[patchID];

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
                  - (nCurrent & pSigma)
                  + impK*(n & pGradDD)
                )*rImpK
            )
        );
    }
}

  
void pRveNonLinGeomTotalLagSolid::updateAvgSigma(const bool debug)
{
    boundBox box(mesh().points());
        
    scalar totalVolume = box.span().x()*box.span().y()*box.span().z();

        // forAll(mesh().boundary(), patchI)
        // {
        //     if (isA<cyclicFvPatch>(mesh().boundary()[patchI]))
        //     {
        //         totalVolume +=
        //             gSum
        //             (
        //                 mesh().Cf().boundaryField()[patchI]
        //               & mesh().Sf().boundaryField()[patchI]
        //             )/2;
        //     }
        // }
        
    scalar Vs = gSum(mesh().V().field());

    if (debug)
    {
        Info << "Vs = " << Vs << ", totalVolume = " << totalVolume << endl;
    }

    // Calculate second Piola-Kirchhoff stress field by transformation
    // of Cauchy stress field
    tensorField PI =
        J_.internalField()
       *(Finv_.internalField() & sigma().internalField());
    
    // Calculate average second Piola-Kirchhoff stress
    avgP_.value() = gSum(PI*mesh().V().field())/totalVolume;

    // Transform average second Piola-kirchhoff stress to Cauchy stress
    tensor curAvgF = (avgRelF_.value() & avgF_.value());
    avgSigma_.value() =
        symm(avgP_.value() & curAvgF.T())/det(curAvgF);
    
    // avgSigma_.value() =
    //     gSum(sigma().internalField()*mesh().V().field())
    //    /(totalVolume + SMALL);

    if (debug)
    {
        Info << "Avg Cauchy stress: " << avgSigma_.value() << endl;
    }
}

  
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
