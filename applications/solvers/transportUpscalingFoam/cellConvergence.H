{
    currCellIter++;

    volTensorField gradX(fvc::grad(X));
    
    if(deltaDRes)
    {
        //- Calculate effective diffusion coefficient
        dimensionedTensor Deff
        (
            "Deff",
            fvc::domainIntegrate
            (
                D*beta*(tensor::I + fvc::grad(X))
                +  ((beta*Ueff-Ustar)*X) 

            ) / vol
        );


        //- Calculate alternative diffusion coefficient
        dimensionedTensor DeffII
        (
            "DeffII",
            fvc::domainIntegrate
            (
                beta*D
                * (tensor::I + gradX)
                &
                T(tensor::I + gradX)
            )  / vol
        );



            XRes = mag(tr(Deff).value() -tr(DeffII).value())
                /(mag(tr(Deff).value())+small);
    }
    else
    {
        XRes = fvc::domainIntegrate(mag(X -X.prevIter())).value()
            /(fvc::domainIntegrate(mag(X)).value()+small);
    }

    cellConverged =
        (
            XRes < corrTol
        )
        ||
        (
            currCellIter > maxCellProbIter
        );

    Info<<"Residual : " << XRes << " Converged : " << cellConverged
        <<endl;

}
