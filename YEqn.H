tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.schemes().div("div(phi,Yi_h)")
    )
);


Info << "Recalculating chemistry" << endl;
reaction->correct();

forAll(selectedScalarFields, i)
{
    const word& fieldName = selectedScalarFields[i];
    Info << "Writing RR." + fieldName << endl;
    volScalarField::Internal& RRfield(mesh.lookupObjectRef<volScalarField::Internal>("RR." + fieldName));
    RRfield.write();
}


if (writeQdot)
{
    Info << "Writing Qdot field" << endl;
    volScalarField Qdot(reaction->Qdot());
    Qdot.write();
}

if (writeFlameSpeed)
{
    volScalarField flameSpeed
    (
        IOobject
        (
            "flameSpeed",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimVelocity
    );
    volScalarField diffusionTerm
    (
        IOobject
        (
            "diffusionTerm",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimMass/dimVolume/dimTime
    );
    volScalarField sourceTerm
    (
        IOobject
        (
            "sourceTerm",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimMass/dimVolume/dimTime
    );

    thermo.correct();

    forAll(Y, i)
    {
        scalarField source(Y[i].size());
        scalarField diffusion(Y[i].size());

        if (Y[i].name() == "O2")
        {
            volScalarField& Yi = Y[i];
            volScalarField magGradY(mag(fvc::grad(Yi)));
            magGradY.max(1e-12);
            //diffusionTerm = fvc::laplacian(thermophysicalTransport->DEff(Yi), Yi);
            diffusionTerm = fvc::laplacian(thermo.alphahe(), Yi);

            // get unburnt background density
            List<scalar> startDistance(Pstream::nProcs(), 0.);
            List<scalar> rhoStartOfEachCore(Pstream::nProcs());
            startDistance[Pstream::myProcNo()] = mag(mesh.C()[0] - vector::zero);
            rhoStartOfEachCore[Pstream::myProcNo()] = rho[0];
            Pstream::gatherList(startDistance);
            Pstream::gatherList(rhoStartOfEachCore);
            Pstream::scatterList(startDistance);
            Pstream::scatterList(rhoStartOfEachCore);
            scalar unburntRho = rhoStartOfEachCore[findMax(startDistance)];

            forAll(Yi, ele)
            {
                diffusion[ele] = diffusionTerm[ele];
                sourceTerm[ele] = (-reaction->R(Yi)->source()[ele])/mesh.V()[ele];
                flameSpeed[ele] = (sourceTerm[ele] + diffusion[ele])/(magGradY[ele])/unburntRho;
            }
            flameSpeed.write();
            diffusionTerm.write();
            sourceTerm.write();
            magGradY.write();
            
        }

    }

    composition.normalise();
}


