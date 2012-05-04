{
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf(fvc::interpolate(rAU));

    U = rAU*UEqn.H();
    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rAU, rho, U, phi)
    );

    adjustPhi(phiU, U, p);

    phi = phiU +
    (
	    fvc::interpolate(rho)*(g & mesh.Sf())
      + (fvc::interpolate(fc) & mesh.Sf())
      - fvc::snGrad(pc) * mesh.magSf()

    )*rAUf;

    //relax
    if (!pimple.finalIter()) {
        phi = 0.7 * phi_old + 0.3 * phi;
    } else {
	    phi_old = phi;
    }

    for (int nonOrth=0; nonOrth<=pimple.nNonOrthCorr(); nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUf, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));

        pEqn.solve
        (
            mesh.solver(p.select(pimple.finalInnerIter(corr, nonOrth)))
        );

        if (nonOrth == pimple.nNonOrthCorr())
        {
            phi -= pEqn.flux();
        }
    }

    U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();

    #include "continuityErrs.H"

    /*p == p_rgh + rho*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }*/
}
