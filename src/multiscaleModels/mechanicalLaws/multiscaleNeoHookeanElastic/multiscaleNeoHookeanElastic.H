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

Class
    multiscaleNeoHookeanElastic

Description
    Neo-Hookean elasticity, similar to model proposed in Simo & Hughes (1998)
    in Eqn 9.2.6:

        tau = J*p*I + mu*dev(J^(-2.0/3.0)*(F & F.T()))

    where
    J        Jacobian of the deformation gradient
    p        hydrostatic pressure
    I        indentity tensor
    mu       shear modulus
    F        deformation gradient

    If thermal stresses are active then 3.0*K*alpha*deltaT*I is added.

    Reference
    Simo & Hughes, Computational Inelasticity, 1998, Springer.


SourceFiles
    multiscaleNeoHookeanElastic.C

Author
    Zeljko Tukovi, FSB Zagreb. All rights reserved.
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef multiscaleNeoHookeanElastic_H
#define multiscaleNeoHookeanElastic_H

#include "mechanicalLaw.H"
#include "solidModel.H"
#include "pRveNonLinGeomTotalLagTotalDispSolid.H"
#include "pRveNonLinGeomUpdatedLagSolid.H"
#include "pRveNonLinGeomTotalLagSolid.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearElastic Declaration
\*---------------------------------------------------------------------------*/

class multiscaleNeoHookeanElastic
:
    public mechanicalLaw
{
    // Private data
   
        //- Apply rve at the boundary patches
        bool applyBoundaryRve_;
    
        //- Rve solvers for cells and boundary faces
        mutable PtrList<solidModels::pRveNonLinGeomTotalLagSolid> pRveSolvers_;

        // Initial density
        const dimensionedScalar rho_;

        // Young's modulus
        const dimensionedScalar E_;

        // Poisson's ratio
        const dimensionedScalar nu_;

        // Shear modulus
        const dimensionedScalar mu_;

        // Bulk modulus
        const dimensionedScalar K_;
    
        // Young's modulus - patch
        const dimensionedScalar pE_;

        // Shear modulus - patch
        const dimensionedScalar pmu_;

        // Bulk modulus - patch
        const dimensionedScalar pK_;

        // Total deformation gradient volField
        volTensorField* FPtr_;

        // Total deformation gradient surfaceField
        surfaceTensorField* FfPtr_;

    // Private Member Functions
    
        //- Make the F field
        void makeF();

        //- Return a reference to the F field
        volTensorField& F();

        //- Make the Ff field
        void makeFf();

        //- Return a reference to the Ff field
        surfaceTensorField& Ff();

        //- Disallow default bitwise copy construct
        multiscaleNeoHookeanElastic(const multiscaleNeoHookeanElastic&);

        //- Disallow default bitwise assignment
        void operator=(const multiscaleNeoHookeanElastic&);

public:

    //- Runtime type information
    TypeName("multiscaleNeoHookeanElastic");

    // Static data members


    // Constructors

        //- Construct from dictionary
        multiscaleNeoHookeanElastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );


    // Destructor

        virtual ~multiscaleNeoHookeanElastic();


    // Member Functions

        //- Return density
        virtual tmp<volScalarField> rho() const;

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Calculate the stress
        virtual void correct(volSymmTensorField& sigma);

        //- Calculate the stress
        virtual void correct(surfaceSymmTensorField& sigma);
  
        //- Update the yield stress: called at end of time-step
        virtual void updateTotalFields();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
