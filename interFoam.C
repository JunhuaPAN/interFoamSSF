/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 Anton Kidess
     \\/     M anipulation  | A.Kidess@tudelft.nl
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

Application
    interFoamSSF

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach.

    The momentum and other fluid properties are of the "mixture" and a single
    momentum equation is solved.

    interFoam is extended with the SSF formulation to reduce spurious currents.
    See Raeini, Blunt, Bijeljic 2012
    http://dx.doi.org/10.1016/j.jcp.2012.04.011

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties/interfaceProperties.H"
#include "twoPhaseMixture.H"
#include "interpolationTable.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        twoPhaseProperties.correct();

        rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;

        // set up subcycling
        Foam::TimeState pts = runTime.subCycle(2);
        Foam::TimeState ss_pts = ++runTime;

        if (pimple.nOuterCorr() < 2)
        {
            Info<< "ABORT: Outer correctors < 2" << endl;
            return -1;
        }

        // --- Pressure-velocity PIMPLE corrector loop
        for (pimple.start(); pimple.loop(); pimple++)
        {
            runTime.TimeState::operator=(ss_pts);
            surfaceScalarField rhoPhiSum(0.0*rhoPhi);
         	if (pimple.corr() == 0) //should move this outside of the PIMPLE loop
            {
                Info << "Subcycle 1" << endl;
		        //update interface location for t = t_n-1 + dt/2
                #include "alphaEqn.H"
                rhoPhiSum += (0.5 * rhoPhi);
                //advance to next half time step
                ss_pts = ++runTime;
                continue;
            } else {
                Info << "Subcycle 2" << endl;
                //update interface location for t_n
                #include "alphaEqn.H"
                //need to accumulate rhoPhi
                rhoPhi = 0.5 * rhoPhi + 0.5 * rhoPhiSum;
            }
            //restore timestep and time for other fields
            ss_pts = runTime;
            runTime.TimeState::operator=(pts);

            //update properties
            twoPhaseProperties.correct();
            rho == alpha1*rho1 + (scalar(1.0) - alpha1)*rho2;
            //update curvature
            interface.correct();
            
            //surface force
            //Cpc should be 0.5 for dynamic problems, > 0.9 for static problems
            scalar Cpc (readScalar
                (
                    alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("Cpc")
                )
            );
            //sharpen interface function
            volScalarField alpha_pc = 1.0/(1.0-Cpc) * 
                (min( max(alpha1,Cpc/2.0), (1.0-Cpc/2.0) ) - Cpc/2.0);
            surfaceScalarField deltasf = fvc::snGrad(alpha_pc);

            surfaceScalarField fcf = interface.sigma()*interface.Kf()*deltasf;
            // relax capillary force
            if (!pimple.finalIter()) {
                fcf = 0.7 * fcf_old + 0.3 * fcf;
            } else {
                fcf_old = fcf;
            }
            // reconstruct for plotting
            fc = fvc::average(fcf*mesh.Sf()/mesh.magSf());
            
            // solve capillary pressure
            fvScalarMatrix pcEqn
            (
                fvm::laplacian(pc) == fvc::div(fcf*mesh.magSf())
            );
            pcEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));
            pcEqn.solve();

            #include "UEqn.H"
            // --- PISO loop
            for (int corr=0; corr<pimple.nCorr(); corr++)
            {
                #include "pEqn.H"
            }

        } // pimple loop

        runTime.endSubCycle();

		Info << "max(U): " << max(U) << endl;
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
