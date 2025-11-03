# Key Formulas

All notation follows the reference paper.

- For a composition $\mu = (\mu_1,\dots,\mu_n)$ we set
  $$
  k_i(\mu) = \#\{j < i : \mu_j > \mu_i\} + \#\{j > i : \mu_j \ge \mu_i\}.
  $$
  The interpolation evaluation point is
  $$
  \widetilde{\mu} = \bigl(q^{\mu_1} t^{-k_1(\mu)}, \dots, q^{\mu_n} t^{-k_n(\mu)}\bigr).
  $$

- The nonsymmetric interpolation Macdonald polynomial $E^*_\mu$ is uniquely determined by
  $$
  [x^\mu]\,E^*_\mu = 1
  \quad\text{and}\quad
  E^*_\mu(\widetilde{\nu}) = 0 \text{ for all } \nu \ne \mu \text{ with } |\nu| \le |\mu|.
  $$

- The type-$A$ Hecke generators act on polynomials in $x_1,\dots,x_n$ via
  $$
  T_i f = t\,f - \frac{t x_i - x_{i+1}}{x_i - x_{i+1}}\,(f - s_i f),
  $$
  where $s_i$ swaps $x_i$ and $x_{i+1}$.

- Given a partition $\lambda$ and a permutation $\sigma_\mu$ sending $\lambda$ to $\mu$, the interpolation ASEP polynomial is
  $$
  f^*_\mu = T_{\sigma_\mu}\,E^*_\lambda.
  $$

- The interpolation Macdonald polynomial is recovered from the signed multiline queue weights as
  $$
  P^*_\lambda(\mathbf{x};q,t) = \sum_{\mu \in S_n(\lambda)} f^*_\mu(\mathbf{x};q,t).
  $$


- The hook normalisation satisfies
  $$
  \text{hook}_\lambda = \prod_{(i,j)\in\lambda} \left(1 - q^{\lambda'_j - i} t^{\lambda_i - j + 1}\right),
  $$
  so that $J^*_\lambda = \text{hook}_\lambda P^*_\lambda$ and $\text{hook}_\lambda f^*_\mu$ are the integral forms.
