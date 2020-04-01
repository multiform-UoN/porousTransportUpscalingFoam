//- Create controls for the power iterations

// Aitken controls
scalar aitkenInitialRes(0);
scalar aitkenInitialResPsi(0);
scalar aitkenInitialResPsiAdj(0);
scalar aitkenTheta(0);
scalar aitkenThetaPsi(0);
scalar aitkenThetaPsiAdj(0);

//scalar aitkenDeltaLambda(0);
//scalar aitkenDeltaLambdaAdj(0);

//- Residual on the eigenvalue
scalar lambdaRes(1.);
scalar psiRes(1.);
scalar psiAdjRes(1.);
scalar XRes(1.);
label currPowerIter(0);
label currCellIter(0);
bool powerConverged(false);
bool cellConverged(false);
simpleControl pwrctrl(mesh,"PowerControl");

const dictionary& pwrDict(pwrctrl.dict());

//- Tolerance for the eigenvalue
scalar eigenValTol
(
    readScalar(pwrDict.lookup("eigenvalueTol"))
);

//- Tolerance for the eigenfunction
scalar eigenFunTol
(
    readScalar(pwrDict.lookup("eigenfunctionTol"))
);

//- Tolerance for the first order corrector
scalar corrTol
(
    readScalar(pwrDict.lookup("cellProblemTol"))
);

//- Max iterations
label maxPowerIter
(
    readLabel(pwrDict.lookup("maxPowerIterations"))
);

label maxCellProbIter
(
    readLabel(pwrDict.lookup("maxCellProblemIterations"))
);

Switch AitkenRelaxX(
    pwrDict.lookupOrDefault<Switch>("AitkenAccelerationCorr",false)
);

if(AitkenRelaxX)
{
	Info<<"\nUsing Aitken acceleration scheme for corrector iterations"
        << endl;
}

Switch AitkenRelaxPsi(
    pwrDict.lookupOrDefault<Switch>("AitkenAccelerationPower",false)
);

if(AitkenRelaxPsi)
{
	Info<<"\nUsing Aitken acceleration scheme for power iterations"
        << endl;
}

Switch implicitRayleigh
(
    pwrDict.lookupOrDefault<Switch>("implicitRayleigh",false)
);

if(implicitRayleigh)
{
    Info <<"\nUsing implicit Rayleigh scheme" << endl;
}

Switch computeResolvent(
    pwrDict.lookupOrDefault<Switch>("computeResolvent",false)
);
