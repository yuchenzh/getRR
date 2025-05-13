/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2022 OpenFOAM Foundation
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

Application
    reactingFoam

Description
    Transient solver for turbulent flow of compressible reacting fluids with
    optional mesh motion and mesh topology changes.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fluidReactionThermo.H"
#include "combustionModel.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidReactionThermophysicalTransportModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

	timeSelector::addOptions();
    argList::addOption
    (
        "fields",
        "list",
        "specify a list of species to be processed. Eg, '(NH3 NO)' - "
        "regular expressions not currently supported"
    );

    argList::addBoolOption
    (
        "Qdot",
        "specify whether to write the heat release rate field. "
    );

    argList::addBoolOption
    (
        "flameSpeed",
        "specify whether to write the flame speed field. "
    );
    //#include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    List<word> selectedScalarFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedScalarFields;
    }

    bool writeQdot = args.optionFound("Qdot");
    bool writeFlameSpeed = args.optionFound("flameSpeed");

    Info << "Selected fields: " << selectedScalarFields << endl;
    Info << "Write Qdot: " << writeQdot << endl;
    Info << "Write flame speed: " << writeFlameSpeed << endl;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI],timeI);
        #include "createMesh.H"
        #include "createDyMControls.H"
        #include "initContinuityErrs.H"
        #include "createFields.H"
        //#include "createFieldRefs.H"
        #include "createRhoUfIfPresent.H"
        Info << "Time = " << runTime.timeName() << nl << endl;

        #include "YEqn.H"
        
        
    }





    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
