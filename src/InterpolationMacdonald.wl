(* ::Package:: *)

BeginPackage["InterpolationMacdonald`"];

TildeVector::usage =
  "TildeVector[mu, q, t] gives the Knop–Sahi interpolation point \\!\\(\\*SubscriptBox[\\(\\~\\), \\(\\)]\\)\\(\\mu\\) for the composition mu.";

NonsymmetricInterpolationMacdonald::usage =
  "NonsymmetricInterpolationMacdonald[mu, vars, {q, t}] returns the nonsymmetric interpolation Macdonald polynomial E^*_mu in the variables vars.";

HeckeOperator::usage =
  "HeckeOperator[poly, i, vars, {q, t}] applies the type-A Hecke operator T_i to the polynomial poly.";

InterpolationASEP::usage =
  "InterpolationASEP[mu, vars, {q, t}] evaluates the interpolation ASEP polynomial f^*_mu via Hecke operators.";

InterpolationASEPFamily::usage =
  "InterpolationASEPFamily[lambda, vars, {q, t}] computes f^*_mu for every permutation mu of the partition lambda.";

InterpolationMacdonaldPolynomial::usage =
  "InterpolationMacdonaldPolynomial[lambda, vars, {q, t}] returns the interpolation Macdonald polynomial P^*_lambda.";

InterpolationElementaryStar::usage =
  "InterpolationElementaryStar[k, vars, t] returns the inhomogeneous elementary symmetric function e_k^*(vars; t) described in Section 8 (q = 1 factorisation).";

InterpolationASEPPartialSum::usage =
  "InterpolationASEPPartialSum[lambda, subset, vars, t] sums f^*_mu(vars; q = 1, t) over all permutations of lambda whose support is the given subset.";

InterpolationASEPPartialProduct::usage =
  "InterpolationASEPPartialProduct[lambda, subset, vars, t] evaluates the right-hand side of the q = 1 factorisation in Theorem 8.1 for the given subset.";

InterpolationMacdonaldQOne::usage =
  "InterpolationMacdonaldQOne[lambda, vars, t] computes P^*_lambda(vars; q = 1, t) using the factorised product of e^*_k terms.";

NonsymmetricInterpolationMacdonald::arg =
  "Variable list length `2` must match the composition length `1`.";
HeckeOperator::arg =
  "Hecke generator index `1` is incompatible with `2` variables.";
InterpolationASEP::orbit =
  "Composition `1` is not in the permutation orbit of `2`.";

Begin["`Private`"];

ClearAll[
  kVector, tildePower, compositionsExact, compositionsUpToDegree,
  monomialFromExponent, solveInterpolationSystem, permutationForOrbit,
  permutationReducedWord, applyHeckeWord, normalizeParams, normalizeT,
  subsetFactor, supportPermutations, partitionConjugate];

kVector[mu_List] := Module[{n = Length[mu]}, Table[
    Count[mu[[;; i - 1]], _?(# > mu[[i]] &)] +
    Count[mu[[i + 1 ;;]], _?(# >= mu[[i]] &)],
    {i, 1, n}
  ]];

tildePower[mu_List, q_, t_] := Module[{k = kVector[mu]},
   q^mu*t^(-k)
];

TildeVector[mu_List, q_: q, t_: t] := tildePower[mu, q, t];

(* -- compositions utilities ------------------------------------------------ *)

normalizeParams[Automatic] := {Symbol["q"], Symbol["t"]};
normalizeParams[{a_, b_}] := {a, b};
normalizeParams[_] := {Symbol["q"], Symbol["t"]};

normalizeT[Automatic] := Symbol["t"];
normalizeT[t_] := t;

compositionsExact[0, 0] := {{}};
compositionsExact[0, _Integer?Positive] := {};
compositionsExact[_, s_Integer?Negative] := {};
compositionsExact[n_Integer?Positive, s_Integer?NonNegative] :=
  compositionsExact[n, s] = Module[{result = {}},
    Do[
      result = Join[result,
        Map[Prepend[#, k] &, compositionsExact[n - 1, s - k]]
      ],
      {k, 0, s}
    ];
    result
  ];

compositionsUpToDegree[n_Integer?Positive, d_Integer?NonNegative] :=
  compositionsUpToDegree[n, d] =
    Flatten[Table[compositionsExact[n, s], {s, 0, d}], 1];

monomialFromExponent[vars_List, exp_List] :=
  Times @@ MapThread[#1^#2 &, {vars, exp}];

solveInterpolationSystem[mu_List, vars_List, {q_, t_}] := Module[
  {n = Length[mu], d = Total[mu], allExponents, unknownExponents, unknownVars,
   evaluationCompositions, matrix, rhs, tilde, coord, coeffs, monomials},
  
  allExponents = compositionsUpToDegree[n, d];
  evaluationCompositions = DeleteCases[allExponents, mu];
  unknownExponents = DeleteCases[allExponents, mu];
  If[unknownExponents === {},
    Return[monomialFromExponent[vars, mu]]
  ];
  unknownVars = Array[c, Length[unknownExponents]];
  
  matrix = Table[
    tilde = TildeVector[nu, q, t];
    monomialFromExponent[tilde, alpha],
    {nu, evaluationCompositions}, {alpha, unknownExponents}
  ];
  
  rhs = -Table[
     tilde = TildeVector[nu, q, t];
     monomialFromExponent[tilde, mu],
     {nu, evaluationCompositions}
   ];
  
  coeffs = LinearSolve[matrix, rhs];
  
  monomials = Map[monomialFromExponent[vars, #] &, unknownExponents];
  Expand[
    monomialFromExponent[vars, mu] + monomials.coeffs
  ]
];

(* -- public E^* function --------------------------------------------------- *)

NonsymmetricInterpolationMacdonald[mu_List, vars_: Automatic, params_: Automatic] :=
  NonsymmetricInterpolationMacdonald[mu, vars, params] =
    Module[{varList, n = Length[mu], qSym, tSym},
      {qSym, tSym} = normalizeParams[params];
      varList = Replace[vars, Automatic :> Array[x, n]];
      If[Length[varList] =!= n,
        Message[NonsymmetricInterpolationMacdonald::arg, n, Length[varList]];
        Return[$Failed];
      ];
      solveInterpolationSystem[mu, varList, {qSym, tSym}]
    ];

(* -- Hecke operators ------------------------------------------------------- *)

HeckeOperator::usageDetail =
  "The polynomial is assumed to use the same ordering of the variables passed via vars.";

HeckeOperator[poly_, i_Integer?Positive, vars_: Automatic, params_: Automatic] :=
  Module[{varList, n, tSym, xi, xip1, swapRules, sPoly},
    varList = Replace[vars, Automatic :> Array[x, Max[i + 1, 2]]];
    n = Length[varList];
    If[i >= n,
      Message[HeckeOperator::arg, i, n];
      Return[$Failed];
    ];
    tSym = normalizeParams[params][[2]];
    xi = varList[[i]];
    xip1 = varList[[i + 1]];
    swapRules = Thread[{xi, xip1} -> {xip1, xi}];
    sPoly = poly /. swapRules;
    Expand[tSym*poly - ((tSym*xi - xip1)/(xi - xip1))*(poly - sPoly)]
  ];

(* -- permutation helpers --------------------------------------------------- *)

permutationForOrbit[lambda_List, mu_List] := Module[
  {n = Length[lambda], positions = <||>, sigmaInv, sigma},
  If[Sort[lambda] =!= Sort[mu],
    Message[InterpolationASEP::orbit, mu, lambda];
    Return[$Failed];
  ];
  Do[
    positions[lambda[[i]]] = Append[Lookup[positions, lambda[[i]], {}], i],
    {i, n}
  ];
  sigmaInv = Table[
    Module[{val = mu[[i]], available},
      available = Lookup[positions, mu[[i]], {}];
      If[available === {},
        Message[InterpolationASEP::orbit, mu, lambda];
        Return[$Failed, Module];
      ];
      positions[val] = Rest[available];
      First[available]
    ],
    {i, n}
  ];
  sigmaInv
];

permutationReducedWord[perm_List] := Module[
  {list = perm, moves = {}, n = Length[perm], pos},
  Do[
    pos = First@First@Position[list, i];
    While[pos < i,
      AppendTo[moves, pos];
      {list[[pos]], list[[pos + 1]]} = {list[[pos + 1]], list[[pos]]};
      pos++;
    ],
    {i, n, 1, -1}
  ];
  Reverse[moves]
];

applyHeckeWord[poly_, word_List, vars_, params_] :=
  Fold[HeckeOperator[#1, #2, vars, params] &, poly, word];

subsetFactor[subset_List, vars_List, tVal_] := Module[
  {sorted = Sort[subset], n = Length[vars]},
  Times @@ Table[
     vars[[sorted[[j]]]] - tVal^(sorted[[j]] - n - (j - 1)),
     {j, Length[sorted]}
   ]
];

supportPermutations[lambda_List, subset_List] := Module[
  {n = Length[lambda], sortedSubset = Sort[subset], positivesCount,
   positives, zeroCount, complementLen, perms},
  If[!SubsetQ[Range[n], sortedSubset], Return[{}]];
  positivesCount = LengthWhile[lambda, # > 0 &];
  positives = Take[lambda, positivesCount];
  zeroCount = n - positivesCount;
  complementLen = n - Length[sortedSubset];
  If[Length[sortedSubset] =!= positivesCount || complementLen =!= zeroCount,
    Return[{}]
  ];
  perms = DeleteDuplicates[Permutations[positives]];
  Table[
    With[{perm = perms[[idx]]},
      Module[{mu = ConstantArray[0, n]},
        Do[
          mu[[sortedSubset[[j]]]] = perm[[j]],
          {j, Length[sortedSubset]}
        ];
        mu
      ]
    ],
    {idx, Length[perms]}
  ]
];

partitionConjugate[lambda_List] := Module[
  {max = If[lambda === {}, 0, Max[lambda]]},
  Table[Count[lambda, _?(# >= r &)], {r, 1, max}]
];

(* -- q = 1 factorisation utilities ---------------------------------------- *)

InterpolationElementaryStar[k_Integer?NonNegative, vars_List, tSym_: Automatic] :=
  Module[{n = Length[vars], tVal = normalizeT[tSym]},
    Which[
      k < 0 || k > n, 0,
      True,
      Total[
        subsetFactor[#, vars, tVal] & /@ Subsets[Range[n], {k}]
      ]
    ]
  ];

InterpolationASEPPartialSum[lambda_List, subset_List, vars_: Automatic, tSym_: Automatic] :=
  Module[{varList = Replace[vars, Automatic :> Array[x, Length[lambda]]],
    tVal = normalizeT[tSym], mus, sums, n, positivesCount, zeroCount,
    complementLen},
    n = Length[varList];
    positivesCount = LengthWhile[lambda, # > 0 &];
    zeroCount = n - positivesCount;
    complementLen = n - Length[subset];
    If[Length[subset] =!= positivesCount || complementLen =!= zeroCount,
      Return[0]
    ];
    mus = supportPermutations[lambda, subset];
    sums = InterpolationASEP[#, varList, {1, tVal}] & /@ mus;
    Expand[Total[sums]]
  ];

InterpolationASEPPartialProduct[lambda_List, subset_List, vars_: Automatic, tSym_: Automatic] :=
  Module[{varList = Replace[vars, Automatic :> Array[x, Length[lambda]]],
    tVal = normalizeT[tSym], lambdaPrime, firstFactor, secondFactor, n,
    positivesCount, zeroCount, complementLen},
    n = Length[varList];
    positivesCount = LengthWhile[lambda, # > 0 &];
    zeroCount = n - positivesCount;
    complementLen = n - Length[subset];
    If[Length[subset] =!= positivesCount || complementLen =!= zeroCount,
      Return[0]
    ];
    lambdaPrime = partitionConjugate[lambda];
    firstFactor = subsetFactor[subset, varList, tVal];
    secondFactor = If[Length[lambdaPrime] <= 1,
      1,
      Times @@ (InterpolationElementaryStar[#, varList, tVal] & /@ Rest[lambdaPrime])
    ];
    Expand[firstFactor*secondFactor]
  ];

InterpolationMacdonaldQOne[lambda_List, vars_: Automatic, tSym_: Automatic] :=
  Module[{varList = Replace[vars, Automatic :> Array[x, Length[lambda]]],
    tVal = normalizeT[tSym], lambdaPrime, factors},
    lambdaPrime = partitionConjugate[lambda];
    factors = InterpolationElementaryStar[#, varList, tVal] & /@ lambdaPrime;
    Expand[If[factors === {}, 1, Times @@ factors]]
  ];

(* -- interpolation ASEP and Macdonald polynomials ------------------------- *)

InterpolationASEP[mu_List, vars_: Automatic, params_: Automatic] :=
  Module[{lambda, permutation, word, base, normParams},
    lambda = Sort[mu, Greater];
    normParams = normalizeParams[params];
    base = NonsymmetricInterpolationMacdonald[lambda, vars, normParams];
    permutation = permutationForOrbit[lambda, mu];
    If[permutation === $Failed, Return[$Failed]];
    word = permutationReducedWord[permutation];
    applyHeckeWord[base, word, Replace[vars, Automatic :> Array[x, Length[mu]]], normParams] // Expand
  ];

InterpolationASEPFamily[lambda_List, vars_: Automatic, params_: Automatic] :=
  Module[{orbit = DeleteDuplicates@Permutations[lambda], normParams},
    normParams = normalizeParams[params];
    AssociationThread[orbit, InterpolationASEP[#, vars, normParams] & /@ orbit]
  ];

InterpolationMacdonaldPolynomial[lambda_List, vars_: Automatic, params_: Automatic] :=
  Module[{orbit = DeleteDuplicates@Permutations[lambda], normParams},
    normParams = normalizeParams[params];
    Total[InterpolationASEP[#, vars, normParams] & /@ orbit] // Expand
  ];

End[];

EndPackage[];
