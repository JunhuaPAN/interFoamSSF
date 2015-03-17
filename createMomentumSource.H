Info<< "\nRead chi" << endl;

IOdictionary chiProperties
(
    IOobject
    (
        "chi",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const scalar chi(readScalar(chiProperties.lookup("chi")));

vector vSource( chi * g.value() * rho2.value() );

volVectorField fSourceMom
(
    IOobject
    (
        "fSourceMom",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("fSourceMom0", dimMass/(dimLength*dimLength*dimTime*dimTime), vector(0,0,0))
);

forAll(mesh.cells(),cellI)
{
	for(label iDir=0; iDir<3; iDir++) fSourceMom[cellI][iDir] = vSource[iDir];
}

Info << "Momemtum source " << chi << " g * rho2 (liqiud) = " << vSource << endl;

surfaceScalarField fSourceMomf ( fvc::interpolate(fSourceMom) & mesh.Sf() );
