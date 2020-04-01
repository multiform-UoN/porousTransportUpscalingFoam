if(AitkenRelaxPsi)
{
	// aitken acceleration
	if(currPowerIter == 0)
	{
		Info << "Updating Aitken initial residuals" << endl;
		aitkenInitialResPsi = max(mag(psi.internalField())).value();
		aitkenInitialResPsiAdj = max(mag(psiAdj.internalField())).value();
	}

	aitkenDeltaPsi.storePrevIter();
	aitkenDeltaPsiAdj.storePrevIter();

	// update delta
	aitkenDeltaPsi = (psi - psi.prevIter()) / aitkenInitialResPsi;
	aitkenDeltaPsiAdj = (psiAdj - psiAdj.prevIter()) / aitkenInitialResPsiAdj;

	// update relaxation factor
	if(currPowerIter == 0)
	{
	    aitkenThetaPsi = 0.1;
		aitkenThetaPsiAdj = 0.1;
	}
	else
	{
	    scalarField b
	    (
	  	   aitkenDeltaPsi.internalField()
	     - aitkenDeltaPsi.prevIter().internalField()
	    );

	    scalar sumMagB = gSum(magSqr(b));
	    if(sumMagB < small)
	    {
	        Info << "Warning: Aitken under-relaxation: denominator less then small"
	             << endl;
	        sumMagB += small;
	    }

	    aitkenThetaPsi = -aitkenThetaPsi*
	        gSum(aitkenDeltaPsi.prevIter().internalField() * b)/sumMagB;

	    Info << "Aitken relaxation factor psi : "<< aitkenThetaPsi << endl;

		b = (
			   aitkenDeltaPsiAdj.internalField()
			  - aitkenDeltaPsiAdj.prevIter().internalField()
	    );

	    sumMagB = gSum(magSqr(b));
	    if(sumMagB < small)
	    {
	        Info << "Warning: Aitken under-relaxation: denominator less then small"
	             << endl;
	        sumMagB += small;
	    }

	    aitkenThetaPsiAdj = -aitkenThetaPsiAdj*
	        gSum(aitkenDeltaPsiAdj.prevIter().internalField() * b)/sumMagB;

	    Info << "Aitken relaxation factor psiAdj : "<< aitkenThetaPsiAdj << endl;

	}

	// correction to the latest psi and psiAdj
	psi += aitkenThetaPsi*aitkenDeltaPsi*aitkenInitialResPsi;
	psiAdj += aitkenThetaPsiAdj*aitkenDeltaPsiAdj*aitkenInitialResPsiAdj;
}
// if(AitkenRelax)
// {
// 	// aitken acceleration
// 	if(currPowerIter == 0)
// 	{
// 	    aitkenInitialResPsi =mag(lambda.value() - lambdaOld);
// 		aitkenInitialResPsiAdj =mag(lambdaAdj.value() - lambdaAdjOld);
// 	}
//
// 	scalar adlOld(aitkenDeltaLambda);
// 	scalar adlaOld(aitkenDeltaLambdaAdj);
//
// 	// update delta
// 	aitkenDeltaLambda = (lambda.value() - lambdaOld) / aitkenInitialResPsi;
// 	aitkenDeltaLambdaAdj = (lambdaAdj.value() - lambdaAdjOld) / aitkenInitialResPsiAdj;
//
// 	// update relaxation factor
// 	if(currPowerIter == 0)
// 	{
// 	    aitkenThetaPsi = 0.1;
// 		aitkenThetaPsiAdj = 0.1;
// 	}
// 	else
// 	{
// 	    scalar b
// 	    (
// 	  	   aitkenDeltaLambda
// 	     - adlOld
// 	    );
//
// 	    scalar sumMagB = magSqr(b);
// 	    if(sumMagB < small)
// 	    {
// 	        Info << "Warning: Aitken under-relaxation: denominator less then small"
// 	             << endl;
// 	        sumMagB += small;
// 	    }
//
// 	    aitkenThetaPsi = -aitkenThetaPsi*adlOld*b/sumMagB;
//
// 	    Info << "Aitken relaxation factor lambda : "<< aitkenThetaPsi << endl;
//
// 		b = (
// 			aitkenDeltaLambdaAdj
//  	     - adlaOld
// 	    );
//
// 	    sumMagB = magSqr(b);
// 	    if(sumMagB < small)
// 	    {
// 	        Info << "Warning: Aitken under-relaxation: denominator less then small"
// 	             << endl;
// 	        sumMagB += small;
// 	    }
//
// 	    aitkenThetaPsiAdj = -aitkenThetaPsiAdj*adlaOld*b/sumMagB;
//
// 	    Info << "Aitken relaxation factor lambdaAdj : "<< aitkenThetaPsiAdj << endl;
//
// 	}
//
// 	// correction to the latest psi and psiAdj
// 	lambda.value() += aitkenThetaPsi*aitkenDeltaLambda*aitkenInitialResPsi;
// 	lambdaAdj.value() += aitkenThetaPsiAdj*aitkenDeltaLambdaAdj*aitkenInitialResPsiAdj;
// }
