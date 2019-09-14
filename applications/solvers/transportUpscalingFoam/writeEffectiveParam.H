Info<<"**********************************************************"<<endl
    <<"                    EFFECTIVE PARAMETERS                "<<endl
    <<"\n   Dispersion tensor Symm : "
    << DeffII.value()<<endl
    <<"\n   Dispersion tensor Symm (trace) : "
    << tr(DeffII.value())<<endl
    <<"\n   Dispersion tensor : "
    << Deff.value()<<endl
    <<"\n   Dispersion tensor  (trace) : "
    << tr(Deff.value())<<endl
    <<"\n   Dispersion tensor (scaled) : "
    << (Deff/D).value()<<endl
    <<"\n   Lambda " << lambda.value() << endl
    <<"\n   First order corrector coefficient : "
    << CI.value() << endl
    <<"\n   Effective velocity : "
    << Ueff.value() << endl;

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
      << (DeffII).value().yy() << ", "<< Ueff.value().x();
