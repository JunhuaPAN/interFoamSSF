/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "interfaceProperties.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvc.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::interfaceProperties::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfaceProperties::correctContactAngle
(
    surfaceVectorField::GeometricBoundaryField& nHatb,
    surfaceVectorField::GeometricBoundaryField& gradAlphaf
) const
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::GeometricBoundaryField& abf = alpha1_.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleFvPatchScalarField& acap =
                const_cast<alphaContactAngleFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );

            fvsPatchVectorField& nHatp = nHatb[patchi];
            const scalarField theta
            (
                convertToRad*acap.theta(U_.boundaryField()[patchi], nHatp)
            );

            const vectorField nf
            (
                boundary[patchi].nf()
            );

            // Reset nHatp to correspond to the contact angle

            const scalarField a12(nHatp & nf);
            const scalarField b1(cos(theta));

            scalarField b2(nHatp.size());
            forAll(b2, facei)
            {
                b2[facei] = cos(acos(a12[facei]) - theta[facei]);
            }

            const scalarField det(1.0 - a12*a12);

            scalarField a((b1 - a12*b2)/det);
            scalarField b((b2 - a12*b1)/det);

            nHatp = a*nf + b*nHatp;
            nHatp /= (mag(nHatp) + deltaN_.value());

            acap.gradient() = (nf & nHatp)*mag(gradAlphaf[patchi]);
            acap.evaluate();
        }
    }
}


void Foam::interfaceProperties::calculateK()
{
/*    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();*/

    scalar CSK = 0.5;
    volScalarField alpha1s_ = CSK * (fvc::average(fvc::interpolate(alpha1_))) + (1.0 - CSK) * alpha1_;
    volScalarField alpha2s_ = CSK * fvc::average(fvc::interpolate(alpha1s_)) + (1.0 - CSK) * alpha1s_;
    volScalarField alpha3s_ = CSK * fvc::average(fvc::interpolate(alpha2s_)) + (1.0 - CSK) * alpha2s_;

    volVectorField gradAlpha = fvc::grad(alpha3s_);
    volVectorField ns(gradAlpha/(Foam::mag(gradAlpha) + deltaN_));

    //clip alpha1 to avoid sqrt(<0)
    volScalarField alpha1c_ = Foam::max(1.0, Foam::min(alpha1_, 0.0));

    // eq. 14 in the reference paper does not have the minus sign,
    // but that's inconsistent with interFoam's formulation
    K_ = -fvc::div(ns);

    volScalarField w = Foam::sqrt(alpha1_*(1.0 - alpha1_) + 1e-6);
    //volScalarField factor = 2.0*Foam::sqrt(alpha1_*(1.0 - alpha1_)+1e-6) * pos(alpha1_-1e-6) * pos(0.9999999-alpha1_);
    volScalarField factor = 2.0 * Foam::sqrt(alpha1c_*(1.0 - alpha1c_));
    //volScalarField factor = 2.0 * w; // alternative to clipping alpha1 - use w

    volScalarField Ks_star_ = fvc::average(fvc::interpolate(K_*w))/fvc::average(fvc::interpolate(w));
    volScalarField Ks1_ = factor * K_ + (1.0 - factor) * Ks_star_;

    Ks_star_ = fvc::average(fvc::interpolate(Ks1_*w))/fvc::average(fvc::interpolate(w));
    volScalarField Ks2_ = factor * K_ + (1.0 - factor) * Ks_star_;

    Kf_ = fvc::interpolate(w*Ks2_)/fvc::interpolate(w);


    // Cell gradient of alpha
    /*const volVectorField gradAlpha(fvc::grad(alpha1_));
    // Interpolated face-gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));*/

    /*gradAlphaf -=
        (mesh.Sf()/mesh.magSf())
       *(fvc::snGrad(alpha1_) - (mesh.Sf() & gradAlphaf)/mesh.magSf());*/

    // Face unit interface normal
/*    surfaceVectorField nHatfv(gradAlphaf/(mag(gradAlphaf) + deltaN_));
    correctContactAngle(nHatfv.boundaryField(), gradAlphaf.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatfv & Sf;

    // Simple expression for curvature
    K_ = -fvc::div(nHatf_);
*/
    // Complex expression for curvature.
    // Correction is formally zero but numerically non-zero.
    /*
    volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
    forAll(nHat.boundaryField(), patchi)
    {
        nHat.boundaryField()[patchi] = nHatfv.boundaryField()[patchi];
    }

    K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
    */
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProperties::interfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    cAlpha_
    (
        readScalar
        (
            alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("cAlpha")
        )
    ),
    sigma_(dict.lookup("sigma")),

    deltaN_
    (
        "deltaN",
        1e-8/pow(average(alpha1.mesh().V()), 1.0/3.0)
    ),

    alpha1_(alpha1),
    U_(U),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("nHatf", dimArea, 0.0)
    ),

    K_
    (
        IOobject
        (
            "K",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("K", dimless/dimLength, 0.0)
    ),

    Kf_
    (
        IOobject
        (
            "Kf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("Kf", dimless/dimLength, 0.0)
    )

{
    calculateK();
}


// ************************************************************************* //
