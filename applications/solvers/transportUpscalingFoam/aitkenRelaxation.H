if(AitkenRelaxX)
{
	// aitken acceleration
	if(currCellIter == 0)
	{
	    aitkenInitialRes = max(mag(X.internalField())).value();
	}

	aitkenDelta.storePrevIter();

	// update delta
	aitkenDelta = (X - X.prevIter()) / aitkenInitialRes;

	// update relaxation factor
	if(currCellIter == 0)
	{
	    aitkenTheta = 0.1;
	}
	else
	{
	    vectorField b
	    (
	  	   aitkenDelta.internalField()
	     - aitkenDelta.prevIter().internalField()
	    );

	    scalar sumMagB = gSum(magSqr(b));
	    if(sumMagB < small)
	    {
	        Info << "Warning: Aitken under-relaxation: denominator less then small"
	             << endl;
	        sumMagB += small;
	    }

	    aitkenTheta = -aitkenTheta*
	        gSum(aitkenDelta.prevIter().internalField() & b)/sumMagB;

	    Info << "Aitken relaxation factor : "<< aitkenTheta << endl;
	}

	// correction to the latest DU
	X += aitkenTheta*aitkenDelta*aitkenInitialRes;
}
