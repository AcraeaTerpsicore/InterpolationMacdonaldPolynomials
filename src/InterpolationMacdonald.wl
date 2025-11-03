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

HookProduct::usage =
  "HookProduct[lambda, {q, t}] computes \\!\\(\\hook_\\lambda\\) for the partition lambda.";

IntegralInterpolationMacdonald::usage =
  "IntegralInterpolationMacdonald[lambda, vars, {q, t}] returns the integral form polynomial J^*_lambda.";

IntegralInterpolationASEP::usage =
  "IntegralInterpolationASEP[mu, vars, {q, t}] multiplies f^*_mu by \\!\\(\\hook_\\lambda\\) where lambda is the partition sorted from mu.";

SignedQueueTableauQ::usage =
  "SignedQueueTableauQ[lambda, tableau] checks the defining constraints of a signed queue tableau determined by lambda.";

SignedQueueTableauStatistics::usage =
  "SignedQueueTableauStatistics[lambda, tableau, {q, t}] returns the statistics used in the tableau weight formula.";

SignedQueueTableauPolynomial::usage =
  "SignedQueueTableauPolynomial[lambda, tableau, vars, {q, t}] evaluates \\!\\(\\mathrm{wt}(\\phi) x^\\phi\\) for a tableau.";

SignedQueueTableauWeight::usage =
  "SignedQueueTableauWeight[lambda, tableau, vars, {q, t}] returns an association describing the tableau weight factor, monomial, and combined contribution.";

SignedQueueTableauType::usage =
  "SignedQueueTableauType[lambda, tableau, n] determines the composition type defined in Section 7.";

SignedMultilineQueueQ::usage =
  "SignedMultilineQueueQ[queue] recognises the strand-based representation of signed multiline queues.";

SignedMultilineQueueFromTableau::usage =
  "SignedMultilineQueueFromTableau[lambda, tableau] constructs the signed multiline queue (strand data) obtained from the tableau via Tab.";

SignedMultilineQueueToTableau::usage =
  "SignedMultilineQueueToTableau[queue] reconstructs the tableau data associated with a signed multiline queue.";

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
  subsetFactor, supportPermutations, partitionConjugate, hookProduct,
  tableauData, tableauLambda, tableauRows, tableauColumns, tableauCellExistsQ,
  tableauValue, tableauRowType, tableauLevel, tableauDownRow, tableauUpRow,
  tableauRestrictedQ, tableauUnrestrictedQ, tableauTopEntryValues,
  signedQueueTableauAttacks, signedQueueTableauMaj, signedQueueTableauCoinv,
  signedQueueTableauEmp, signedQueueTableauNegative, signedQueueTableauStats,
  SignedMultilineQueueFromTableau, SignedMultilineQueueToTableau,
  SignedMultilineQueueQ];

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

tableauLambda[lambda_List] := Replace[lambda, {} -> {0}];
tableauRows[lambda_List] := 2*Max[tableauLambda[lambda]];
tableauColumns[lambda_List] := Length[lambda];

tableauData[lambda_List, tableau_] := Module[{data},
  data = Which[
    MatchQ[tableau, _Association?(KeyExistsQ[#, "Data"] &)], tableau["Data"],
    True, tableau
  ];
  If[!ListQ[data], Return[$Failed]];
  If[Length[data] =!= tableauRows[lambda], Return[$Failed]];
  If[!And @@ (ListQ[#] && Length[#] === tableauColumns[lambda] & /@ data),
    Return[$Failed]
  ];
  data
];

tableauCellExistsQ[lambda_List, row_Integer, col_Integer] :=
  With[{height = If[col <= Length[lambda], lambda[[col]], 0]},
    height > 0 && row <= 2*height
  ];

tableauValue[data_List, row_Integer, col_Integer] := Module[{cols = Length[data[[1]]]},
  If[row < 1 || row > Length[data] || col < 1 || col > cols,
    Null,
    data[[row, col]]
  ]
];

tableauRowType[row_Integer] := If[OddQ[row], "Classic", "Primed"];
tableauLevel[row_Integer] := If[OddQ[row], (row + 1)/2, row/2];
tableauDownRow[row_Integer] := If[row <= 1, None, row - 1];
tableauUpRow[row_Integer, maxRow_Integer] := If[row >= maxRow, None, row + 1];

tableauRestrictedQ[lambda_List, data_List, row_Integer, col_Integer] := Module[
  {down = tableauDownRow[row], value = tableauValue[data, row, col], downValue},
  Which[
    value === Null, True,
    down === None, True,
    !tableauCellExistsQ[lambda, down, col], True,
    (downValue = tableauValue[data, down, col]) === Null, True,
    True, Abs[downValue] === Abs[value]
  ]
];

tableauUnrestrictedQ[lambda_List, data_List, row_Integer, col_Integer] :=
  Not[tableauRestrictedQ[lambda, data, row, col]];

tableauTopEntryValues[lambda_List, data_List] := Module[{cols = tableauColumns[lambda]},
  Table[
    With[{height = lambda[[col]], row = 2*lambda[[col]]},
      If[height > 0,
        {height, col, Abs[tableauValue[data, row, col]], tableauValue[data, row, col]},
        Nothing
      ]
    ],
    {col, cols}
  ]
];

signedQueueTableauAttacks[lambda_List, data_List, row_Integer, col_Integer] := Module[
  {cols = tableauColumns[lambda], value = tableauValue[data, row, col], down,
   attacks = {}, rowType = tableauRowType[row], level = tableauLevel[row],
   lambdaCol = If[col <= Length[lambda], lambda[[col]], 0]},
  If[value === Null || value === 0, Return[{}]];
  Do[
    If[c =!= col && tableauCellExistsQ[lambda, row, c] && tableauValue[data, row, c] =!= Null,
      AppendTo[attacks, {row, c}]
    ],
    {c, cols}
  ];
  down = tableauDownRow[row];
  If[down =!= None,
    Which[
      value > 0,
      Do[
        If[c =!= col && tableauCellExistsQ[lambda, down, c] && tableauValue[data, down, c] =!= Null,
          If[c <= Length[lambda] && lambdaCol >= lambda[[c]], AppendTo[attacks, {down, c}]]
        ],
        {c, cols}
      ],
      value < 0,
      Do[
        If[c =!= col && tableauCellExistsQ[lambda, down, c] && tableauValue[data, down, c] =!= Null,
          If[c < col && c <= Length[lambda] && lambda[[c]] > lambdaCol, AppendTo[attacks, {down, c}]]
        ],
        {c, cols}
      ]
    ]
  ];
  attacks
];

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

hookProduct[lambda_List, params_] := Module[
  {lam = Sort[Select[lambda, # > 0 &], Greater], lamPrime, qSym, tSym},
  If[lam === {}, Return[1]];
  lamPrime = partitionConjugate[lam];
  {qSym, tSym} = normalizeParams[params];
  Expand[Times @@ Flatten[
      Table[
        1 - qSym^(lamPrime[[j]] - i) tSym^(lam[[i]] - j + 1),
        {i, Length[lam]}, {j, lam[[i]]}
      ]
    ]]
];

HookProduct[lambda_List, params_: Automatic] :=
  hookProduct[lambda, params];

IntegralInterpolationMacdonald[lambda_List, vars_: Automatic, params_: Automatic] :=
  Module[{normParams = normalizeParams[params]},
    Expand[HookProduct[lambda, normParams] *
      InterpolationMacdonaldPolynomial[lambda, vars, normParams]]
  ];

IntegralInterpolationASEP[mu_List, vars_: Automatic, params_: Automatic] :=
  Module[{normParams = normalizeParams[params], lambda = Sort[mu, Greater]},
    Expand[HookProduct[lambda, normParams] *
      InterpolationASEP[mu, vars, normParams]]
  ];

(* -- signed queue tableaux ------------------------------------------------ *)

signedQueueTableauStats[lambda_List, tableau_, params_: Automatic] :=
  Module[{data = tableauData[lambda, tableau], rows = tableauRows[lambda],
    cols = tableauColumns[lambda], qSym, tSym, restricted, legMatrix, armMatrix,
    classicUnrestricted = {}, primedUnrestricted = {}, maj = 0, coinv = 0,
    emp = 0, negative = 0, topEntries, levels, existsQ, rowValuesClassic,
    columnHeights = lambda, nVars, cellVal, downRow, upRow,
    downVal, attacks, valueList, level, lambdaCol},
    If[data === $Failed, Return[$Failed]];
    {qSym, tSym} = normalizeParams[params];
    restricted = Table[True, {rows}, {cols}];
    legMatrix = Table[0, {rows}, {cols}];
    armMatrix = Table[0, {rows}, {cols}];
    existsQ[row_, col_] := tableauCellExistsQ[columnHeights, row, col];
    rowValuesClassic = Table[{}, {rows}];
    nVars = Length[columnHeights];
    topEntries = tableauTopEntryValues[columnHeights, data];
    (* top entries order check *)
    Do[
      With[{sameHeight = Select[topEntries, First[#] == h &]},
        If[Length[sameHeight] > 1,
          If[!And @@ (sameHeight[[#, 3]] > sameHeight[[# + 1, 3]] & /@ Range[Length[sameHeight] - 1]),
            Return[$Failed]
          ]
        ]
      ],
      {h, DeleteDuplicates[columnHeights]}
    ];
    (* cell validations *)
    Do[
      If[existsQ[row, col],
        cellVal = tableauValue[data, row, col];
        If[!IntegerQ[cellVal] || cellVal == 0, Return[$Failed]];
        If[Abs[cellVal] > nVars, Return[$Failed]];
        If[tableauRowType[row] === "Classic" && cellVal <= 0, Return[$Failed]];
        downRow = tableauDownRow[row];
        If[downRow =!= None && existsQ[downRow, col],
          downVal = tableauValue[data, downRow, col];
          If[tableauRowType[row] === "Primed" && downVal =!= Null,
            If[Abs[downVal] < Abs[cellVal], Return[$Failed]]
          ];
          If[cellVal < 0 && downVal =!= Null && Abs[downVal] =!= Abs[cellVal],
            negative++
          ];
        ];
        If[tableauRowType[row] === "Classic",
          rowValuesClassic[[row]] = Append[rowValuesClassic[[row]], cellVal]
        ];
        restricted[[row, col]] = tableauRestrictedQ[columnHeights, data, row, col];
        If[tableauRowType[row] === "Classic",
          legMatrix[[row, col]] = columnHeights[[col]] - tableauLevel[row],
          legMatrix[[row, col]] = columnHeights[[col]] - tableauLevel[row] + 1
        ];
        If[tableauRowType[row] === "Classic" && restricted[[row, col]] === False,
          AppendTo[classicUnrestricted, {row, col}]
        ];
        If[tableauRowType[row] === "Primed" && restricted[[row, col]] === False,
          AppendTo[primedUnrestricted, {row, col}]
        ];
      , (* else cell not in diagram *)
        If[tableauValue[data, row, col] =!= Null, Return[$Failed]]
      ],
      {row, rows}, {col, cols}
    ];
    (* positive in primed rows must appear in classic row *)
    Do[
      If[existsQ[row, col],
        cellVal = tableauValue[data, row, col];
        If[tableauRowType[row] === "Primed" && cellVal > 0,
          downRow = tableauDownRow[row];
          If[downRow === None || !MemberQ[rowValuesClassic[[downRow]], cellVal],
            Return[$Failed]
          ]
        ]
      ],
      {row, rows}, {col, cols}
    ];
    (* attacking constraints *)
    Do[
      If[existsQ[row, col],
        cellVal = tableauValue[data, row, col];
        attacks = signedQueueTableauAttacks[columnHeights, data, row, col];
        If[!And @@ (Abs[cellVal] =!= Abs[tableauValue[data, ##[[1]], ##[[2]]]] & /@ attacks),
          Return[$Failed]
        ]
      ],
      {row, rows}, {col, cols}
    ];
    (* maj statistic *)
    Do[
      If[existsQ[row, col] && tableauRowType[row] === "Classic", 
        cellVal = tableauValue[data, row, col];
        downRow = tableauDownRow[row];
        If[downRow =!= None && existsQ[downRow, col],
          downVal = tableauValue[data, downRow, col];
          If[downVal =!= Null && Abs[downVal] < cellVal,
            maj += legMatrix[[row, col]] + 1
          ]
        ]
      ],
      {row, rows}, {col, cols}
    ];
    (* compute arm matrix for classic cells *)
    Do[
      If[existsQ[row, col] && tableauRowType[row] === "Classic" &&
         Not[restricted[[row, col]]],
        lambdaCol = columnHeights[[col]];
        armMatrix[[row, col]] =
          Total[Boole[columnHeights[[#]] < lambdaCol] & /@
             Range[col + 1, cols]] +
          Total[Boole[columnHeights[[#]] == lambdaCol &&
                Not[restricted[[row, #]]]] & /@ Range[col + 1, cols]];
      ],
      {row, rows}, {col, cols}
    ];
    (* arm for primed cells *)
    Do[
      If[existsQ[row, col] && tableauRowType[row] === "Primed" &&
         Not[restricted[[row, col]]],
        upRow = tableauUpRow[row, rows];
        If[upRow =!= None && existsQ[upRow, col] &&
           Not[restricted[[upRow, col]]],
          armMatrix[[row, col]] = armMatrix[[upRow, col]],
          lambdaCol = columnHeights[[col]];
          level = tableauLevel[row];
          armMatrix[[row, col]] =
            Total[Boole[columnHeights[[#]] >= level] & /@
               Range[col + 1, cols]] +
            Total[Boole[columnHeights[[#]] == lambdaCol &&
                 upRow =!= None && upRow <= rows &&
                 Not[restricted[[upRow, #]]]] & /@ Range[1, col - 1]];
        ];
      ],
      {row, rows}, {col, cols}
    ];
    (* coinversion statistic *)
    Do[
      If[existsQ[row, col],
        cellVal = tableauValue[data, row, col];
        If[cellVal > 0,
          downRow = tableauDownRow[row];
          If[downRow =!= None && existsQ[downRow, col],
            downVal = Abs[tableauValue[data, downRow, col]];
            Do[
              If[existsQ[downRow, c],
                valueList = tableauValue[data, downRow, c];
                If[valueList =!= Null,
                  lambdaCol = columnHeights[[col]];
                  level = columnHeights[[c]];
                  upRow = tableauUpRow[downRow, rows];
                  If[(columnHeights[[c]] < lambdaCol) ||
                     (columnHeights[[c]] == lambdaCol && upRow =!= None &&
                      existsQ[upRow, c] && Not[restricted[[upRow, c]]]),
                    With[{absX = cellVal, absY = Abs[valueList]},
                      If[(absX < absY < downVal) ||
                          (downVal < absX < absY) ||
                          (absY < downVal < absX),
                        coinv++
                      ]
                    ]
                  ]
                ]
              ],
              {c, col + 1, cols}
            ]
          ]
        ]
      ],
      {row, rows}, {col, cols}
    ];
    (* empty statistic *)
    Do[
      If[existsQ[row, col] && tableauRowType[row] === "Primed",
        cellVal = tableauValue[data, row, col];
        downRow = tableauDownRow[row];
        If[downRow =!= None && existsQ[downRow, col],
          downVal = tableauValue[data, downRow, col];
          If[downVal =!= Null,
            With[{a = Abs[cellVal], cVal = downVal, classicVals = rowValuesClassic[[downRow]]},
              If[cVal > a,
                emp += Total[Boole[!MemberQ[classicVals, #]] & /@ Range[a + 1, cVal - 1]]
              ]
            ];
          ];
        ];
      ],
      {row, rows}, {col, cols}
    ];
    <|
      "Data" -> data,
      "Restricted" -> restricted,
      "LegMatrix" -> legMatrix,
      "ArmMatrix" -> armMatrix,
      "ClassicUnrestricted" -> classicUnrestricted,
      "PrimedUnrestricted" -> primedUnrestricted,
      "Maj" -> maj,
      "Coinv" -> coinv,
      "Emp" -> emp,
      "Negative" -> negative
    |>
  ];

SignedQueueTableauQ[lambda_List, tableau_] :=
  signedQueueTableauStats[lambda, tableau] =!= $Failed;

SignedQueueTableauStatistics[lambda_List, tableau_, params_: Automatic] :=
  Module[{stats = signedQueueTableauStats[lambda, tableau, params]},
    If[stats === $Failed, Return[$Failed]];
    KeyTake[stats, {"Maj", "Coinv", "Emp", "Negative", "ClassicUnrestricted",
      "PrimedUnrestricted", "LegMatrix", "ArmMatrix", "Restricted"}]
  ];

SignedQueueTableauPolynomial[lambda_List, tableau_, vars_: Automatic, params_: Automatic] :=
  Module[{weight = SignedQueueTableauWeight[lambda, tableau, vars, params]},
    If[weight === $Failed, Return[$Failed]];
    weight["Polynomial"]
  ];

SignedQueueTableauWeight[lambda_List, tableau_, vars_: Automatic, params_: Automatic] :=
  Module[{stats = signedQueueTableauStats[lambda, tableau, params], qSym, tSym,
    varList, nVars, factorClassic, factorPrimed, monomial, rows, cols, data,
    primedCells},
    If[stats === $Failed, Return[$Failed]];
    {qSym, tSym} = normalizeParams[params];
    data = stats["Data"];
    rows = tableauRows[lambda];
    cols = tableauColumns[lambda];
    varList = Replace[vars, Automatic :> Array[x, cols]];
    nVars = Length[varList];
    factorClassic = Times @@ (
        (1 - tSym)/(1 - qSym^(stats["LegMatrix"][[Sequence @@ #]] + 1) *
            tSym^(stats["ArmMatrix"][[Sequence @@ #]] + 1)) & /@
          stats["ClassicUnrestricted"]);
    factorPrimed = (1 - tSym)^(Length[stats["PrimedUnrestricted"]]);
    primedCells = Flatten[Table[
        If[EvenQ[row] && tableauCellExistsQ[lambda, row, col], {row, col}, Nothing],
        {row, rows}, {col, cols}], 1];
    monomial = Times @@ (Function[{rc},
         Module[{row = rc[[1]], col = rc[[2]], val = data[[rc[[1]], rc[[2]]]], level},
           If[val === Null, 1,
             level = tableauLevel[row];
             If[val > 0,
               varList[[val]],
               -qSym^(level - 1)/tSym^(nVars - 1)
             ]
           ]
         ]
       ] /@ primedCells);
    With[{negative = stats["Negative"], maj = stats["Maj"], coinv = stats["Coinv"],
       emp = stats["Emp"]},
      Module[{factor = (-1)^negative * qSym^maj * tSym^(coinv + emp) *
          factorClassic * factorPrimed},
        <|
          "Factor" -> factor,
          "Monomial" -> monomial,
          "Polynomial" -> Expand[factor*monomial],
          "Statistics" -> KeyTake[stats, {"Maj", "Coinv", "Emp", "Negative"}]
        |>
      ]
    ]
  ];

SignedQueueTableauType[lambda_List, tableau_, n_Integer?Positive] :=
  Module[{data = tableauData[lambda, tableau], cols = tableauColumns[lambda],
    heights = lambda, result},
    If[data === $Failed, Return[$Failed]];
    result = ConstantArray[0, n];
    Do[
      If[tableauCellExistsQ[lambda, 1, col],
        With[{val = data[[1, col]]},
          If[IntegerQ[val] && 1 <= val <= n,
            result[[val]] = heights[[col]]
          ]
        ]
      ],
      {col, cols}
    ];
    result
  ];

SignedMultilineQueueFromTableau::nolambda =
  "Unable to infer a partition from the tableau; please supply lambda explicitly.";

SignedMultilineQueueFromTableau[tableau_] :=
  (Message[SignedMultilineQueueFromTableau::nolambda]; $Failed);

SignedMultilineQueueFromTableau[lambda_List, tableau_] :=
  Module[{data = tableauData[lambda, tableau], rows = tableauRows[lambda],
    cols = tableauColumns[lambda]},
    If[data === $Failed, Return[$Failed]];
    <|
      "Lambda" -> lambda,
      "Strands" ->
        Table[
          With[{entries =
             Cases[Table[{row, data[[row, col]]}, {row, rows}],
               {row_, val_} /; val =!= Null :>
                 <|"Row" -> row, "Value" -> val|>]},
            <|"Column" -> col, "Entries" -> entries|>
          ],
          {col, cols}
        ]
    |>
  ];

SignedMultilineQueueToTableau[queue_Association] :=
  Module[{lambda = Lookup[queue, "Lambda", $Failed], rows, cols, data,
    strands},
    If[lambda === $Failed || !VectorQ[lambda, IntegerQ[#] && # >= 0 &],
      Return[$Failed]
    ];
    rows = tableauRows[lambda];
    cols = tableauColumns[lambda];
    data = Table[ConstantArray[Null, cols], {rows}];
    strands = Lookup[queue, "Strands", {}];
    MapIndexed[
      With[{colIndex = Lookup[#1, "Column", First[#2]],
        entries = Lookup[#1, "Entries", {}]},
        With[{boundedCol =
           If[cols == 0, 1,
             Min[Max[colIndex, 1], cols]]},
          Scan[
            With[{row = Lookup[#, "Row", 0], val = Lookup[#, "Value", Null]},
              If[IntegerQ[row] && 1 <= row <= Max[rows, 1] &&
                 val =!= Null,
                If[rows > 0 && cols > 0,
                  data[[row, boundedCol]] = val
                ]
              ]
            ]&,
            entries
          ]
        ]
      ]&,
      strands
    ];
    data
  ];

SignedMultilineQueueQ[queue_] :=
  Module[{lambda, strands, rows, cols, structuralOK},
    If[!AssociationQ[queue] || !KeyExistsQ[queue, "Lambda"] ||
       !KeyExistsQ[queue, "Strands"],
      Return[False]
    ];
    lambda = queue["Lambda"];
    If[!VectorQ[lambda, IntegerQ[#] && # >= 0 &], Return[False]];
    strands = queue["Strands"];
    If[!ListQ[strands], Return[False]];
    rows = tableauRows[lambda];
    cols = tableauColumns[lambda];
    structuralOK = And @@ MapIndexed[
       With[{strand = #1, idx = #2[[1]]},
         With[{col = Lookup[strand, "Column", idx],
           entries = Lookup[strand, "Entries", {}]},
           IntegerQ[col] &&
           (cols == 0 || (1 <= col <= cols)) &&
           ListQ[entries] &&
           And @@ (AssociationQ[#] &&
               IntegerQ[Lookup[#, "Row", 0]] &&
               (rows == 0 || 1 <= Lookup[#, "Row", 0] <= rows) &&
               IntegerQ[Lookup[#, "Value", 0]] & /@ entries)
         ]
       ]&,
       strands
     ];
    If[!structuralOK, Return[False]];
    SignedQueueTableauQ[lambda, SignedMultilineQueueToTableau[queue]]
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
