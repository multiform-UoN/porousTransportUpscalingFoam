/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019
     \\/     M anipulation  | Matteo Icardi, Federico Municchi
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    specCellFoam

Description
    Compute leading order and first correction of the homogenised advection
    diffusion equation with reactive internal boundaries by solving the
    cell problem with spectral decomposition.

Developers
    Federico Municchi, Nottingham (2019)
    Matteo Icardi,  Nottingham (2019)
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallDist.H"
#include "simpleControl.H"


int main(int argc, char *argv[]) {

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createCellProblemControls.H"

    #include "createFields.H"

    // Since solver contains no time loop it would never execute
    // function objects so do it ourselves
    runTime.functionObjects().start();
    // Manually advance the time index
    runTime++;

    runTime.write();

    //- Power iterations (i.e., spectral problem)
    Info << "Power iterations for the spectral problem" << endl;
    do
    {

        #include "powerIterSettings.H"

        while(pwrctrl.correctNonOrthogonal())
        {

            #include    "psiEqn.H"
            #include    "psiAdjEqn.H"

            if(pwrctrl.finalNonOrthogonalIter())
            {

                #include    "updateEigenvalue.H"
                #include    "normalisePsi.H"
                #include    "normalisePsiPsiAdj.H"
                #include    "aitkenRelaxEigenfunctions.H"
                #include    "powerConvergence.H"
            }
        }

    } while (
        !powerConverged
    );

    runTime.write();

    //- Calculate modified velocity and beta
    #include    "calculateBeta.H"
    #include    "calculateUstar.H"

    //- Solve corrector problem
    Info << "Solving for scalar transport corrector field" << endl;
    do {

        #include    "cellIterSettings.H"

        while(pwrctrl.correctNonOrthogonal())
        {
            #include    "XEqn.H"
        }

        #include "rescaleX.H"
        #include "aitkenRelaxation.H"
        #include "cellConvergence.H"

    } while (
        !cellConverged
    );

    #include    "rescaleX.H"
    #include    "calculateEffectiveParam.H"

    runTime.functionObjects().end();
    runTime.writeAndEnd();

    #include    "writeEffectiveParam.H"

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End" << endl;

    return 0;
}
