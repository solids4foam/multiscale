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

#include "pRveNonLinGeomUpdatedLagSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "bound.H"
#include "cyclicFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(pRveNonLinGeomUpdatedLagSolid, 0);
addToRunTimeSelectionTable
(
    solidModel, pRveNonLinGeomUpdatedLagSolid, dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pRveNonLinGeomUpdatedLagSolid::pRveNonLinGeomUpdatedLagSolid
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
    // (
    //     IOobject
    //     (
    //         "avgRelF",
    //         runTime.timeName(),
    //         mesh(),
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh(),
    //     dimensionedTensor
    //     (
    //         "I",
    //         dimless,
    //         // I
    //         tensor(1, 0, 0, -2.75331e-05, 1.00003, -2.23866e-22, 0, 0, 1)
    //     )
    // ),
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
    relF_
    (
        IOobject
        (
            "relF",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        I + gradDD().T()
    ),
    relFinv_
    (
        IOobject
        (
            "relFinv",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(relF_)
    ),
    relJ_
    (
        IOobject
        (
            "relJ",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(relF_)
    ),
    rho_
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mechanical().rho()
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
    )
{
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


bool pRveNonLinGeomUpdatedLagSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    int iCorr = 0;
    lduSolverPerformance solverPerfDD;
    blockLduMatrix::debug = 0;

    Info<< "Solving the updated Lagrangian form of the momentum equation for DD"
        << endl;

    // Momentum equation loop
    do
    {
        // Store fields for under-relaxation and residual calculation
        DD().storePrevIter();

        // Momentum equation incremental updated Lagrangian form
        fvVectorMatrix DDEqn
        (
            fvm::d2dt2(rho_, DD())
          + fvc::d2dt2(rho_.oldTime(), D().oldTime())
         == fvm::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          - fvc::laplacian(impKf_, DD(), "laplacian(DDD,DD)")
          + fvc::div(relJ_*relFinv_ & sigma(), "div(sigma)")
          + rho_*g()
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

        // Relative deformation gradient
        relF_ = avgRelF_ + gradDD().T();

        // Inverse relative deformation gradient
        relFinv_ = inv(relF_);

        // Total deformation gradient
        F_ = relF_ & F_.oldTime();

        // Relative Jacobian (Jacobian of relative deformation gradient)
        relJ_ = det(relF_);

        // Jacobian of deformation gradient
        J_ = relJ_*J_.oldTime();

        // Update the momentum equation inverse diagonal field
        // This may be used by the mechanical law when calculating the
        // hydrostatic pressure
        const volScalarField DEqnA("DEqnA", DDEqn.A());

        // Calculate the stress using run-time selectable mechanical law
        gradDD() += (avgRelF_.value().T() - I);
        mechanical().correct(sigma());
        gradDD() -= (avgRelF_.value().T() - I);
    }
    while
    (
       !converged
        (
            iCorr,
            solverPerfDD.initialResidual(),
            solverPerfDD.nIterations(),
            DD()
        )
     && ++iCorr < nCorr()
    );

    // Update gradient of total displacement
    gradD() = fvc::grad(D().oldTime() + DD());

    // Total displacement
    D() = D().oldTime() + DD();

    // Update pointDD as it used by FSI procedure
    mechanical().interpolate(DD(), pointDD());

    // Total displacement at points
    pointD() = pointD().oldTime() + pointDD();

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


tmp<vectorField> pRveNonLinGeomUpdatedLagSolid::tractionBoundarySnGrad
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
    const tensorField& pGradDD = gradDD().boundaryField()[patchID];

    // Patch Cauchy stress
    const symmTensorField& pSigma = sigma().boundaryField()[patchID];

    // Patch relative deformation gradient inverse
    const tensorField& relFinv = relFinv_.boundaryField()[patchID];

    // Patch unit normals (updated configuration)
    const vectorField n = patch.nf();

    // Patch unit normals (deformed configuration)
    vectorField nCurrent = relFinv.T() & n;
    nCurrent /= mag(nCurrent);

    // Testing: let us instead calculate the deformed normals by interpolating
    // displacements to the points and calculating the normals on the deformed
    // patch; as this is how we will actually move the mesh, it will be more
    // consistent.
    // This, however, begs the question: is the cell-centred deformation
    // gradient field 'F' consistent with our point displacement field?"
    // i.e. we can calculate the deformed cell volumes two ways (at least):
    //     1. V = J*Vold
    //     2. Move the mesh with pointD and then directly calculate V
    // The answers from 1. and 2. are only approximately equal: this causes a
    // slight inconsistency. The equalavent can be said for the deformed face
    // areas.
    // In Maneeratana, the mesh is never moved, instead method 1. is used for
    // the deformed volumes and areas.

    // standAlonePatch deformedPatch =
    //     standAlonePatch
    //     (
    //         mesh().boundaryMesh()[patchID].localFaces(),
    //         mesh().boundaryMesh()[patchID].localPoints()
    //     );

    // // Calculate the deformed points
    // const pointField deformedPoints =
    //     mechanical().volToPoint().interpolate
    //     (
    //         mesh().boundaryMesh()[patchID],
    //         DD_
    //     )
    //   + mesh().boundaryMesh()[patchID].localPoints();

    // // Move the standAlonePatch points
    // const_cast<pointField&>(deformedPatch.points()) = deformedPoints;

    // // Patch unit normals (deformed configuration)
    // const vectorField& nCurrent = deformedPatch.faceNormals();

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


void pRveNonLinGeomUpdatedLagSolid::updateTotalFields()
{
    // Density
    rho_ = rho_.oldTime()/relJ_;

    // Move the mesh to the deformed configuration
    const vectorField oldPoints = mesh().allPoints();
    moveMesh(oldPoints, DD(), pointDD());

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
    
    solidModel::updateTotalFields();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
