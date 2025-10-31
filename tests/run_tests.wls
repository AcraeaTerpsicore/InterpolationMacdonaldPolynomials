Needs["InterpolationMacdonald`", "src/InterpolationMacdonald.wl"];

vars2 = Array[x, 2];
vars3 = Array[x, 3];

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

test6 = zeroQ[
   InterpolationASEPPartialProduct[{2, 1, 0}, {1, 3}, vars3, t] - 
     Expand[(vars3[[1]] - t^(-2))*(vars3[[3]] - t^(-1))*
       InterpolationElementaryStar[1, vars3, t]]];

test7 = zeroQ[
   InterpolationMacdonaldQOne[{2, 1, 0}, vars3, t] - 
     Expand[InterpolationElementaryStar[2, vars3, t]*
       InterpolationElementaryStar[1, vars3, t]]];

muList = {{2, 0, 1}, {1, 0, 2}};
test8 = zeroQ[
   InterpolationASEPPartialSum[{2, 1, 0}, {1, 3}, vars3, t] - 
     Total[InterpolationASEP[#, vars3, {1, t}] & /@ muList]];

test9 = zeroQ[
   InterpolationElementaryStar[1, vars2, t] -
     ((vars2[[1]] - t^(-1)) + (vars2[[2]] - 1))];

results = {test1, test2, test3, test4, test5, test6, test7, test8, test9};

If[AllTrue[results, TrueQ],
  Print["TESTS PASSED"],
  Print["TESTS FAILED"];
  Print[results]
];
