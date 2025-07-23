# Simulation Models

Simulation designs from three SDR-related papers are summarised below.
For each paper we list

* **Response formulae**
* **Structural directions** (true \$\beta\$ vectors)
* **Predictor / treatment generation** (Case 1 vs Case 2 when applicable)
* **Dimension settings** for \$p=6\$ versus \$p=12\$
* **Returned values** from the accompanying R helpers
* **Reproducibility checklist**

---

## 1 · Ma & Zhu (2013) – Section 4

### Response Distributions

* **Example 1**

  $$
    Y\mid X \sim N\!\bigl(X^{\!\top}\beta,\;1\bigr)
  $$

* **Example 2**

  $$
    Y\mid X \sim
      N\!\Bigl(\sin\!\bigl(2X^{\!\top}\beta\bigr)+2e^{\,2+X^{\!\top}\beta},
                \ \log\!\bigl\{2+(X^{\!\top}\beta)^2\bigr\}\Bigr)
  $$

* **Example 3**

  $$
    Y\mid X \sim
      N\!\bigl(2(X^{\!\top}\beta_1)^2,\ 2e^{X^{\!\top}\beta_2}\bigr)
  $$

### Structural Directions

$$
\beta   =(1.3,-1.3,1,-0.5,0.5,-0.5)^{\!\top},
\qquad
\beta_1 =(1,\tfrac23,\tfrac23,0,-\tfrac13,\tfrac23)^{\!\top},
\qquad
\beta_2 =(0.8,0.8,-0.3,0.3,0,0)^{\!\top}.
$$

### Predictor Generation

1. Draw $X_1,X_2,e_1,e_2\stackrel{\text{iid}}{\sim}N(0,1)$.
2. $X_3 = 0.2X_1 + 0.2(X_2+2)^2 + 0.2e_1$.
3. $X_4 = 0.1 + 0.1(X_1+X_2) + 0.3(X_1+1.5)^2 + 0.2e_2$.
4. $X_5\sim\text{Bernoulli}\bigl(\operatorname{logit}^{-1}(X_1)\bigr)$.
5. $X_6\sim\text{Bernoulli}\bigl(\operatorname{logit}^{-1}(X_2)\bigr)$.

### Returned Values

`X , Y , beta / beta1 / beta2 , eta  (the appropriate linear index †)`

† `eta = X %*% beta` for Example 1, or `X %*% beta1 / beta2` for Example 3.

### Reproducibility

*Generate $X$ via steps 1–5 → compute indices → sample $Y$ using the formula.*

### Reference

Ma, Y. & Zhu, L. (2013) **Efficient Estimation in Sufficient Dimension Reduction**.

---

## 2 · Ma & Zhu (2012) – Section 6

### Response Formulas  (\$\varepsilon\sim N(0,1)\$)

* Model 1

  $$
    Y=\frac{\eta_1}{0.5+(\eta_2+1.5)^2}+0.5\,\varepsilon
  $$

* Model 2

  $$
    Y=\eta_1^{2}+2\lvert\eta_2\rvert+0.1\lvert\eta_2\rvert\,\varepsilon
  $$

* Model 3

  $$
    Y = e^{\eta_1}+2(\eta_2+1)^{2}+\lvert\eta_1\rvert\,\varepsilon
  $$

* Model 4

  $$
    Y = \eta_1^{2}+\eta_2^{2}+0.5\,\varepsilon
  $$

### Structural Directions

$$
\beta_1=\tfrac1{\sqrt6}(1,1,1,1,1,1,0,\dots,0)^{\!\top},\quad
\beta_2=\tfrac1{\sqrt6}(1,-1,1,-1,1,-1,0,\dots,0)^{\!\top}.
$$

### Predictor Generation

* **Case 1 (Mixed / non-Gaussian)**

  *Gaussian core* — $X_1,X_2\sim N(0,1)$ (plus $X_7{:}X_{12}$ when $p=12$
  from MVN with AR(1) $\rho=0.5$.

  *Additions*

  $$
    \begin{aligned}
      X_3 &= \lvert X_1+X_2\rvert + \lvert X_1\rvert\,e_1,\\
      X_4 &= \lvert X_1+X_2\rvert^{2} + \lvert X_2\rvert\,e_2,\\
      X_5 &\sim\text{Bernoulli}\bigl(\operatorname{logit}^{-1}(X_2)\bigr),\\
      X_6 &\sim\text{Bernoulli}\bigl(\Phi(X_2)\bigr),
    \end{aligned}\quad
    e_{1,2}\sim N(0,1).
  $$

* **Case 2 (Multivariate normal)** — $X\sim N_p\bigl(0,\Sigma\bigr)$ with
  $\Sigma_{ij}=0.5^{|i-j|}$.

### Dimension Settings

* $p=6$: $X_1{:}X_6$ active.
* $p=12$: same first six + Gaussian noise $X_7{:}X_{12}$.

### Returned Values

`X , Y , beta1 , beta2 , eta1 , eta2`

### Reproducibility

1. Choose model (1–4), case (1 or 2), and $p=6$ or 12.
2. Generate $X$, compute $\eta_1, \eta_2$, then sample $Y$.

### Reference

Ma, Y. & Zhu, L. (2012) **A Semiparametric Approach to Dimension Reduction**.

---

## 3 · Nabi et al. (2022) – Section 5

### Components

* **Baseline covariates**
  $C=(C\_1,\dots,C\_4)\sim N\_4(0,I)$
  $C\_{\text{sum}}=\sum C\_i$;
  $C\_{\text{sign}}=\sum(-1)^i C\_i$;
  $C\_{\text{alt}}=C\_1-C\_2-C\_3+C\_4$.

* **Indices**
  $\eta_1 = A^{\top}\beta\_1,\quad \eta\_2 = (A^{\odot2})^{!\top}\beta_2$.

* **Outcome**

  $$
    Y = \eta_1 + \eta_2 + C_{\text{sum}}
        + \bigl(\textstyle\sum_j A_j\bigr)C_{\text{sum}} + \varepsilon,
    \qquad \varepsilon\sim N(0,1).
  $$

### Structural Directions

$$
\beta_1=\tfrac1{\sqrt6}(1,1,1,1,1,1,0,\dots,0)^{\!\top},\quad
\beta_2=\tfrac1{\sqrt6}(1,-1,1,-1,1,-1,0,\dots,0)^{\!\top}.
$$

### Treatment Generation

* **Case 1 (Mixed / non-Gaussian)**

  1. *Gaussian block*

     * $p=6$: generate $(A\_1,A\_2)\$, mean \$(C\_{\text{sum}},C\_{\text{sign}})$, AR(1) $\rho=0.5$.
     * $p=12$: generate $(A\_1,A\_2,A\_7{:}A\_{12})$, mean
       $(C\_{\text{sum}},C\_{\text{sign}},C\_1,C\_2,C\_3,-C\_1,-C\_2,-C\_3)$.
  2. *Non-Gaussian*

     $$
       \begin{aligned}
         A_3 &= \lvert A_1+A_2\rvert + \lvert A_1\rvert\,e_1,\\
         A_4 &= \lvert A_1+A_2\rvert^{2} + \lvert A_2\rvert\,e_2,\\
         A_5 &\sim\text{Bernoulli}\bigl(\operatorname{logit}^{-1}(A_2)\bigr),\\
         A_6 &\sim\text{Bernoulli}\bigl(\Phi(A_2)\bigr),
       \end{aligned}\quad e_{1,2}\sim N(0,1).
     $$

* **Case 2 (Multivariate normal)**
  $A \sim N\_p\bigl(\mu(C),\Sigma\bigr)$ with AR(1) covariance ($\rho=0.5$) and
  row-specific means

  $$
    \mu(C)=
    \begin{cases}
      (C_{\text{sum}},C_{\text{sign}},C_{\text{alt}},-C_{\text{alt}},
       C_{\text{sum}}-2C_3,C_{\text{sum}}-2C_1),
      & p=6,\\[6pt]
      (C_{\text{sum}},C_{\text{sign}},C_{\text{alt}},-C_{\text{alt}},
       C_{\text{sum}}-2C_3,C_{\text{sum}}-2C_1,
       C_1,C_2,C_3,-C_1,-C_2,-C_3),
      & p=12.
    \end{cases}
  $$

### Dimension Settings

* $p=6$: active $A_1{:}A_6$.
* $p=12$: same six + Gaussian noise $A_7{:}A_{12}$.

### Returned Values

`Y , A , C , beta1 , beta2 , eta1 , eta2`

### Reproducibility Checklist

1. Draw $C$, compute $C\_{\text{sum}}$, $C\_{\text{sign}}$, $C\_{\text{alt}}$.
2. Generate $A$ (Case 1 or 2, chosen $p$).
3. Compute $\eta_1,\eta_2$.
4. Sample $Y$.

### Reference

Nabi, R., McNutt, T., & Shpitser, I. (2022) **Semiparametric Causal Sufficient
Dimension Reduction of Multidimensional Treatments**.