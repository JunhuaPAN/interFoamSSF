    
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(twoPhaseProperties.muf(), U)
      - (fvc::grad(U) & fvc::grad(twoPhaseProperties.muf()))
    //- fvc::div(muf*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                fvc::interpolate(rho)*(g & mesh.Sf())
              + (
                  fcf
                  - fvc::snGrad(p)
				  - fvc::snGrad(pc)
                ) * mesh.magSf()
	      // Momemtum source term
              - fSourceMomf	
            )
        );
    }

