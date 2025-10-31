# Interpolation Macdonald Toolkit

This repository implements the constructions from `reference_paper/main.tex` in the Wolfram Language.  The code exposes a small package for working with interpolation Macdonald polynomials via nonsymmetric interpolation polynomials and Hecke operators.

## Layout

- `src/InterpolationMacdonald.wl` – main package with helpers for interpolation points, Hecke operators, interpolation ASEP polynomials, and the symmetric interpolation Macdonald polynomials.
- `tests/run_tests.wls` – automated regression script executed with `wolframscript`.
- `FORMULAS.md` – key definitions and weights transcribed from the reference paper.
- `TEST_SUMMARY.md` – latest test results.

## Getting Started

Load the package inside Mathematica or a `wolframscript` session:

```wolfram
Needs["InterpolationMacdonald`", "src/InterpolationMacdonald.wl"];
InterpolationASEP[{1, 0}, Array[x, 2]]
```

All exported functions accept optional parameter lists `{q, t}`.  When omitted, the symbols `q` and `t` from the global context are used automatically.

## Running the Tests

The tests require `wolframscript` (Mathematica 14).  From the project root run:

```bash
"/mnt/d/Software/Wolfram Research/Mathematica/14.0/wolframscript.exe" -file tests/run_tests.wls
```

`TEST_SUMMARY.md` is kept in sync with the most recent run.

## Example

```wolfram
Needs["InterpolationMacdonald`", "src/InterpolationMacdonald.wl"];
vars = Array[x, 2];
f10 = InterpolationASEP[{1, 0}, vars];
P11 = InterpolationMacdonaldPolynomial[{1, 1}, vars];
```

This snippet reproduces the signed multiline queue formulae for the smallest non-trivial compositions and partitions.
