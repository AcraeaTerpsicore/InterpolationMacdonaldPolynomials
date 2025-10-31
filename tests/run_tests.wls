Needs["InterpolationMacdonald`", "src/InterpolationMacdonald.wl"];

vars2 = Array[x, 2];

zeroQ[expr_] := TrueQ[PossibleZeroQ[Simplify[expr]]];

test1 = zeroQ[InterpolationASEP[{0, 0}, vars2] - 1];
test2 = zeroQ[
   InterpolationASEP[{1, 0}, vars2] - 
     NonsymmetricInterpolationMacdonald[{1, 0}, vars2]];

topComponent[poly_, vars_List, deg_Integer] := 
  SeriesCoefficient[poly /. Thread[vars -> vars*s], {s, 0, deg}];

macTop = topComponent[
   InterpolationMacdonaldPolynomial[{1, 1}, vars2], vars2, 2];
fFamily = InterpolationASEPFamily[{1, 1}, vars2];
macTopFromFamily = 
  Total[topComponent[#, vars2, 2] & /@ Values[fFamily]];
test3 = zeroQ[macTop - macTopFromFamily];

tilde00 = TildeVector[{0, 0}, q, t];
test4 = zeroQ[InterpolationASEP[{1, 0}, vars2] /. 
    Thread[vars2 -> tilde00]];

test5 = zeroQ[
   Total[Values[fFamily]] - 
     InterpolationMacdonaldPolynomial[{1, 1}, vars2]];

results = {test1, test2, test3, test4, test5};

If[AllTrue[results, TrueQ],
  Print["TESTS PASSED"],
  Print["TESTS FAILED"];
  Print[results]
];
