Needs["InterpolationMacdonald`", "src/InterpolationMacdonald.wl"];

vars2 = Array[x, 2];
vars3 = Array[x, 3];

zeroQ[expr_] := TrueQ[PossibleZeroQ[Simplify[expr]]];
tDenominatorQ[den_] := den === 1 ||
   Module[{vars = Variables[den], exp},
     vars === {t} && (exp = Exponent[den, t]; exp >= 0 &&
        TrueQ[Simplify[den/t^exp] == 1])
   ];
noDenomQ[expr_] := tDenominatorQ[Denominator[Together[expr]]];

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

test10 = zeroQ[
   HookProduct[{1}, {q, t}] - (1 - t)];

test11 = noDenomQ[IntegralInterpolationMacdonald[{1, 1}, vars2]];

test12 = noDenomQ[IntegralInterpolationASEP[{1, 0}, vars2]];

lambdaSample = {1, 1};
tableauSample = {{2, 1}, {2, -1}};

test13 = TrueQ[SignedQueueTableauQ[lambdaSample, tableauSample]];

test14 = SignedQueueTableauType[lambdaSample, tableauSample, 2] === {1, 1};

allTableaux =
  Flatten[Table[{{a1, a2}, {b1, b2}}, {a1, {1, 2}}, {a2, {1, 2}}, 
     {b1, {-2, -1, 1, 2}}, {b2, {-2, -1, 1, 2}}], 3];
validTableaux = Select[allTableaux, SignedQueueTableauQ[lambdaSample, #] &];
groupedTableaux = GroupBy[validTableaux, SignedQueueTableauType[lambdaSample, #, 2] &];
test15 = FreeQ[SignedQueueTableauPolynomial[lambdaSample, #, vars2] & /@ validTableaux, $Failed];

statsSample = SignedQueueTableauStatistics[lambdaSample, tableauSample];
test16 = statsSample =!= $Failed && statsSample["Maj"] === 0 &&
   statsSample["Coinv"] === 0 && statsSample["Emp"] === 0 && statsSample["Negative"] === 0;

weightSample = SignedQueueTableauWeight[lambdaSample, tableauSample, vars2];
test17 = weightSample =!= $Failed && 
   SignedQueueTableauPolynomial[lambdaSample, tableauSample, vars2] === 
    weightSample["Polynomial"];

queueSample = SignedMultilineQueueFromTableau[lambdaSample, tableauSample];
test18 = SignedMultilineQueueQ[queueSample];

statsData = InterpolationMacdonald`Private`tableauData[lambdaSample, tableauSample];
test19 = SignedMultilineQueueToTableau[queueSample] === statsData;

queueRoundTrip = SignedMultilineQueueFromTableau[lambdaSample, SignedMultilineQueueToTableau[queueSample]];
test20 = SignedMultilineQueueQ[queueRoundTrip] && queueRoundTrip === queueSample;

test21 = zeroQ[
   SignedQueueTableauPolynomial[lambdaSample, SignedMultilineQueueToTableau[queueSample], vars2] - 
    SignedQueueTableauPolynomial[lambdaSample, tableauSample, vars2]];

queuesEnumerated = SignedMultilineQueueEnumerate[{1, 1}, {1, 1}, vars2];
test22 = zeroQ[
   Total[Lookup[#, "Weight"]["Polynomial"] & /@ queuesEnumerated] - 
    InterpolationASEP[{1, 1}, vars2]];

tableauxEnumerated = SignedQueueTableauEnumerate[{1, 1}, {1, 1}, vars2];
test23 = And @@ {
   zeroQ[Total[Lookup[#, "Weight"]["Polynomial"] & /@ tableauxEnumerated] - 
      InterpolationASEP[{1, 1}, vars2]],
   And @@ (KeyExistsQ[Lookup[#, "Statistics", <||>], "Maj"] & /@ tableauxEnumerated)
   };

pasepSubset = {1, 2};
test24 = zeroQ[
   TPushASEPPartialObservable[{1, 1}, pasepSubset, u, t] - 
    InterpolationASEPPartialProduct[{1, 1}, pasepSubset, {u, u}, t]];

test25 = zeroQ[
   TPushMacdonaldObservable[{2, 1}, vars3, t] - 
    InterpolationMacdonaldQOne[{2, 1}, vars3, t]];

results = {test1, test2, test3, test4, test5, test6, test7, test8, test9,
   test10, test11, test12, test13, test14, test15, test16, test17, test18,
   test19, test20, test21, test22, test23, test24, test25};

If[AllTrue[results, TrueQ],
  Print["TESTS PASSED"],
  Print["TESTS FAILED"];
  Print[results]
];
