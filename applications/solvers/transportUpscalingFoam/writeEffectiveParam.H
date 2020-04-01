Info<<"**********************************************************"<<endl
    <<"                    EFFECTIVE PARAMETERS                "<<endl
    <<"\n   Dispersion tensor Symm : "
    << DeffII.value()<<endl
    <<"\n   Dispersion tensor Symm (trace) : "
    << tr(DeffII.value()) - DeffII.value().zz()<<endl
    <<"\n   Dispersion tensor : "
    << Deff.value()<<endl
    <<"\n   Dispersion tensor  (trace) : "
    << tr(Deff.value()) - Deff.value().zz()<<endl
    <<"\n   Dispersion tensor (scaled) : "
    << (Deff/D).value()<<endl
    <<"\n   Lambda : " << lambda.value() << endl
    <<"\n   LambdaAdj : " << lambdaAdj.value() << endl
    <<"\n   First order corrector coefficient : "
    << CI.value() << endl
    <<"\n   Effective velocity : "
    << Ueff.value() << endl
    <<"\n   Ustar contErrors : "
    << fvc::domainIntegrate(fvc::div(phistar)).value() << endl
    << "\n Advection term effective dispersion : " << endl
    << fvc::domainIntegrate(X*fvc::div(Ustar*X))/vol << endl
    << "\n Diffusion term effective dispersion delta : " << endl
    << D*fvc::domainIntegrate(
      X*fvc::div(beta*(tensor::I + fvc::grad(X)))
      + (beta*fvc::grad(X)&(tensor::I + fvc::grad(X)))
    )/vol << endl
    << "\n Another form of DeffII : "
    << fvc::domainIntegrate(
      D*beta*(tensor::I + fvc::grad(X))
      + beta*D*((tensor::I + fvc::grad(X))&fvc::grad(X))
    )/vol << endl
    << "\n Another relative error : "
    << fvc::domainIntegrate(
      (beta*Ueff-Ustar)*X
    - D*(beta*fvc::grad(X)&(tensor::I + fvc::grad(X)))
  )/vol << endl;

  Info<<"\n Norm beta : " << fvc::domainIntegrate(beta)/vol << endl;

    //- Write dictionary for dispersion tensor
    Info<<"\nCreating effectiveTransportProperties dictionary..."<<endl;
    IOdictionary effectiveTransportProperties
    (
        IOobject
        (
            "effectiveTransportProperties",
            runTime.constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );

    effectiveTransportProperties.add
    (
        "Deff",
        Deff
    );

    effectiveTransportProperties.add
    (
        "DeffSymm",
        DeffII
    );

    effectiveTransportProperties.add
    (
        "lambda",
        lambda.value()
    );

    effectiveTransportProperties.add
    (
        "CI",
        CI
    );

    effectiveTransportProperties.add
    (
        "Ueff",
        Ueff
    );

    effectiveTransportProperties.regIOobject::write();

    //Write matlab format
    fileName outputFile(runTime.path()/"results.csv");
    OFstream os(outputFile);

    // os  << "lambda = " << lambda.value() << "" << endl
    //     << "Dx = " << (Deff/D).value().xx() << "" << endl
    //     << "Ueff = " << Ueff.value().x() << "" << endl;
    os<< lambda.value() << ", " << (DeffII).value().xx() << ", "
      << (DeffII).value().yy() << ", "<< Ueff.value().x()<<", "
      << lambdaAdj.value() << ", "
      << (Deff).value().xx() << ", " << (Deff).value().yy(); 
