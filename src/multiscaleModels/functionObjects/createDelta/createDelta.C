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

\*----------------------------------------------------------------------------*/

#include "createDelta.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(createDelta, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        createDelta,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::createDelta::updateDelta()
{
    // Calculate df
    const labelList& own = mesh_.owner();
    const labelList& nei = mesh_.neighbour();
    const vectorField& C = mesh_.C().internalField();
    
    forAll(df_.internalField(), faceI)
    {
        df_.internalField()[faceI] = (C[nei[faceI]] - C[own[faceI]]);
    }
    
    forAll(df_.boundaryField(), patchI)
    {
        forAll(df_.boundaryField()[patchI], faceI)
        {
            df_.boundaryField()[patchI] =
                mesh_.boundary()[patchI].nf()
               /mesh_.deltaCoeffs().boundaryField()[patchI];
        }
    }  

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::createDelta::createDelta
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    runTime_(t),
    mesh_
    (
        runTime_.lookupObject<fvMesh>
        (
            dict.lookupOrDefault<word>("region", "region0")
        )
    ),
    df_
    (
        IOobject
        (
            "df",
            t.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("0", dimLength, vector::zero)
    )
{
    Info<< "Creating " << this->name()
        << " function object" << endl;
    updateDelta();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::createDelta::start()
{
    return updateDelta();
}

#if FOAMEXTEND > 40
bool Foam::createDelta::execute(const bool forceWrite)
#else
bool Foam::createDelta::execute()
#endif
{
    return updateDelta();
}

bool Foam::createDelta::read(const dictionary& dict)
{
    return true;
}

// ************************************************************************* //
