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
  permutationReducedWord, applyHeckeWord, normalizeParams];

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
  sigma = ConstantArray[0, n];
  Do[
    sigma[[sigmaInv[[i]]]] = i,
    {i, n}
  ];
  sigma
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
