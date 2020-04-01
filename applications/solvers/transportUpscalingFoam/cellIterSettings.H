Info<<endl<<"Cell corrector iteration " << currCellIter <<"\n";
//- Disable screen output from linear solver
//  (Otherwise it looks like the solver never converges)
lduMatrix::debug = 0;

Switch deltaDRes
(
  pwrDict.lookupOrDefault<Switch>("residualsOnEffectiveDispersion",false)
);

X.storePrevIter();
