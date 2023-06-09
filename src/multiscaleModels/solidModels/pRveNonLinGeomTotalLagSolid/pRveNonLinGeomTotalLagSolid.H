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

Class
    pRveNonLinGeomTotalLagSolid

Description
    Mathematical model for finite strains and rotations, where a total
    Lagrangian and incremental displacement approach is employed.

    The stress is calculated by the run-time selectable mechanical law.

    More details in:
    P. Cardiff, Z. Tukovic, P. De Jaeger, M. Clancy and A. Ivankovic. A
    Lagrangian cell-centred finite volume method for metal forming simulation,
    doi=10.1002/nme.5345.

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.
    Philip Cardiff, UCD.  All rights reserved.

SourceFiles
    pRveNonLinGeomTotalLagSolid.C

\*---------------------------------------------------------------------------*/

#ifndef pRveNonLinGeomTotalLagSolid_H
#define pRveNonLinGeomTotalLagSolid_H

#include "solidModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class pRveNonLinGeomTotalLagSolid Declaration
\*---------------------------------------------------------------------------*/

class pRveNonLinGeomTotalLagSolid
:
    public solidModel
{
    // Private data

        //- Average total deformation gradient
        dimensionedTensor avgRelF_;

        //- Average total deformation gradient
        dimensionedTensor avgF_;
    
        //- Average first Piola-Kirchhoff stress
        dimensionedTensor avgP_;
    
        //- Average Cauchy stress
        dimensionedSymmTensor avgSigma_;

        //- Total deformation gradient
        volTensorField F_;

        //- Inverse of the total deformation gradient
        volTensorField Finv_;

        //- Jacobian of total deformation gradient
        volScalarField J_;

        //- Implicit stiffness; coefficient of the Laplacian term
        //  The value of this term only affects convergence and not the answer
        //  In the case of Hookean elasticity, (2*mu + lambda) is the optimum
        const volScalarField impK_;

        //- Implicit stiffness surface field
        surfaceScalarField impKf_;

        //- Reciprocal of impK
        //  We divide by impK every time tractionBoundarySnGrad is called so we
        //  will store the reciprocal for effficiency
        const volScalarField rImpK_;

        //- Cell-to-cell vector
        mutable surfaceVectorField df_;

        //- Time index
        label timeIndex_;
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        pRveNonLinGeomTotalLagSolid(const pRveNonLinGeomTotalLagSolid&);

        //- Disallow default bitwise assignment
        void operator=(const pRveNonLinGeomTotalLagSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::TOTAL_LAGRANGIAN;
        }


public:

    //- Runtime type information
    TypeName("pRveNonLinGeomTotalLag");

    // Constructors

        //- Construct from components
        pRveNonLinGeomTotalLagSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );


    // Destructor

        virtual ~pRveNonLinGeomTotalLagSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                // This model solves for DD
                return DD();
            }

            //- This solid model solves for DD instead of D
            virtual bool incremental() const
            {
                return true;
            }

            // Return current average deformation gradient
            dimensionedTensor& avgRelF()
            {
                if (timeIndex_ != runTime().timeIndex())
                {
                    timeIndex_ = runTime().timeIndex();
                    
                    // Accumulate previous relative deformation gradient
                    avgF_.value() = (avgRelF_.value() & avgF_.value());
                }
                
                return avgRelF_;
            }
  
            // Return current average Cauchy stress
            dimensionedSymmTensor& avgSigma()
            {
                return avgSigma_;
            }

            // Return current average Cauchy stress
            const dimensionedSymmTensor& avgSigma() const
            {
                return avgSigma_;
            }
  
            // Update average rve stress
            void updateAvgSigma(const bool debug = false);
  
        // Edit

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
