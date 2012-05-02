    
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
			fc +
            fvc::reconstruct
            (
                fvc::interpolate(rho)*(g & mesh.Sf())
              + (
                  - fvc::snGrad(p)
				  - fvc::snGrad(pc)
                ) * mesh.magSf()
            )
            /*fvc::reconstruct
            (
                //fvc::interpolate(rho*(scalar(1) - beta*(T-Tref1) ))*(g & mesh.Sf()) +
               (
                  fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1) * factor1 //normal 
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
              + (fvc::interpolate(Fst) & mesh.Sf())
            )*pos(T-Ts1)*/
        );
    }

