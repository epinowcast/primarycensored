# Analytic solutions for censored delay distributions

## 1 Introduction

### 1.1 What are we going to do in this vignette

In this vignette, we’ll derive the analytic solutions for the primary
censored delay distributions for a range of commonly used distributions.

1.  Examine the case of exponentially tilted primary event times
2.  Analyse the uniform primary event time case (\\r=0\\)
3.  Explore general partial expectation and its role in solving equation
    [(3.1)](#eq:unifprim)
4.  Provide specific analytic solutions for distributions where the
    partial expectation can be reduced to an analytic expression

## 2 Analytic solutions for exponentially tilted primary event times

In epidemiological analysis, it is common for primary events to occur at
exponentially increasing or decreasing rates, for example, incidence of
new infections in an epidemic. In this case, the distribution of the
primary event time within its censoring window is biased by the
exponential growth or decay^(\[[1](#ref-Park2024)\]) (i.e. for
exponential growth, the event time will more likely to be closer to the
end of the censoring window and vice versa for decay). If we assume a
reference uniform distribution within a primary censoring window \\\[k,
k + w_P)\\ then the distribution of the primary event time within the
censoring window is the [exponential
tilted](https://en.wikipedia.org/wiki/Exponential_tilting) uniform
distribution:

\\ f_P(t) \propto \exp(r t) \mathbb{1}\_{\[k, k + w_P\]}(t). \\

In this case, the distribution function for \\C_P\\, that is the length
of time left in the primary censor window *after* the primary event
time, is given by:

\\ F\_{C_P}(p; r) = { 1 - \exp(-r p) \over 1 - \exp(-r w_P)}, \qquad p
\in \[0, w_P\]. \tag{2.1} \\ Note that taking the limit \\\lim_r
\rightarrow 0\\ in equation [(2.1)](#eq:fcp) gives the uniform
distribution function \\F\_{C_P}(p, 0) = p / w_P\\.

In the following, it is convenient to use a (linear) difference operator
defined as:

\\ \Delta\_{w}f(t) = f(t + w) - f(t). \\

## 3 Uniform primary event time (\\r=0\\)

Applying a uniform primary event time distribution to equation 3.1 from
[“Why it
works”](https://primarycensored.epinowcast.org/dev/articles/why-it-works.md)
gives:

\\ Q\_{S\_+}(t) = Q_T(t + w_P) + { 1 \over w_P} \int_0^{w_P} f_T(t+p) p~
dp. \\

This is analytically solvable whenever the upper distribution function
of \\T\\ is known and the mean of \\T\\ is analytically solvable from
its integral definition.

In each case considered below it is easier to change the integration
variable:

\\ \begin{aligned} Q\_{S\_+}(t) &= Q_T(t + w_P) + { 1 \over w_P}
\int_t^{t+w_P} f_T(z) (z-t)~ dz \\ &= Q_T(t + w_P) + { 1 \over w_P}
\Big\[ \int_t^{t+w_P} f_T(z) z~ dz - t \Delta\_{w_P}F_T(t) \Big\].
\end{aligned} \tag{3.1} \\

### 3.1 General partial expectation

Note that for any distribution with an analytically available
distribution function \\F_T\\ equation [(3.1)](#eq:unifprim) can be
solved so long as the *partial expectation*

\\ \int_t^{t+w_P} f_T(z) z~ dz \tag{3.2} \\

can be reduced to an analytic expression.

The insight here is that this will be possible for any distribution
where the average of the distribution can be calculated analytically,
which includes commonly used non-negative distributions such as the
Gamma, Log-Normal and Weibull distributions.

### 3.2 General discrete censored delay distribution

First, we note that equation 3.3 from [“Why it
works”](https://primarycensored.epinowcast.org/dev/articles/why-it-works.md)
can be written using the difference operator: \\f_n = -\Delta_1
Q\_{S\_+}(n-1)\\. We can insert this expression into equation
[(3.1)](#eq:unifprim) to give the discrete censored delay distribution
for uniformly distributed primary event times:

\\ \begin{aligned} f_n &= \Delta_1\Big\[(n-1) \Delta_1F_T(n-1)\Big\] -
\Delta_1Q_T(n) - \Delta_1\Big\[ \int\_{n-1}^n f_T(z) z ~dz \Big\] \\ &=
(n+1)F_T(n+1) + (n-1)F_T(n-1) - 2nF_T(n) - \Delta_1\Big\[ \int\_{n-1}^n
f_T(z) z ~dz \Big\]. \end{aligned} \tag{3.3} \\

We now consider some specific delay distributions.

### 3.3 Gamma distributed delay times

The Gamma distribution has the density function:

\\ f_T(z;k, \theta) = {1 \over \Gamma(k) \theta^k} z^{k-1}
\exp(-z/\theta). \\ Where \\\Gamma\\ is the Gamma function.

The Gamma distribution has the distribution function:

\\ \begin{aligned} F_T(z;k, \theta) &= {\gamma(k, z/\theta) \over
\Gamma(k)}, \qquad z\geq 0,\\ F_T(z;k, \theta) &= 0, \qquad z \< 0.
\end{aligned} \\ Where \\\gamma\\ is the lower incomplete gamma
function.

#### 3.3.1 Gamma partial expectation

We know that the full expectation of the Gamma distribution is
\\\mathbb{E}\[T\] = k\theta\\, which can be calculated as a standard
integral. Doing the same integral for the partial expectation gives:

\\ \begin{aligned} \int_t^{t+w_P} f_T(z) z~ dz &= {1 \over \Gamma(k)
\theta^k} \int_t^{t+w_P} \mathbb{1}(z \geq 0) z z^{k-1}
\exp(-z/\theta)~dz \\ &= {\Gamma(k+1) \theta^{k+1} \over \Gamma(k)
\theta^k} {1 \over \Gamma(k+1) \theta^{k+1}} \int_t^{t+w_P} \mathbb{1}(z
\geq 0) z^{k} \exp(-z/\theta)~dz\\ &= k\theta \Delta\_{w_P} F_T(t; k +
1, \theta). \end{aligned} \tag{3.4} \\

#### 3.3.2 Survival function of \\S\_{+}\\ for Gamma distribution

By substituting equation [(3.4)](#eq:gammapartexp) into equation
[(3.3)](#eq:disccensunifprim) we can solve for the survival function of
\\S\_+\\ in terms of analytically available functions:

\\ Q\_{S\_+}(t; k, \theta) = Q_T(t + w_P; k, \theta) + { 1 \over w_P}
\big\[ k \theta \Delta\_{w_P}F_T(t; k+1, \theta) - t \Delta\_{w_P}F_T(t;
k, \theta) \big\]. \tag{3.5} \\

#### 3.3.3 Gamma discrete censored delay distribution

By substituting [(3.5)](#eq:survgammaunifprim) into
[(3.3)](#eq:disccensunifprim) we get the discrete censored delay
distribution in terms of analytically available functions: \\
\begin{aligned} f_n &= (n+1) F_T(n+1; k, \theta) + (n-1) F_T(n-1; k,
\theta) - 2 n F_T(n; k, \theta) - k \theta \Delta_1^{(2)}F_T(n-1; k+1,
\theta)\\ &= (n+1) F_T(n+1; k, \theta) + (n-1) F_T(n-1; k, \theta) - 2 n
F_T(n; k, \theta) \\ &+ k \theta \Big( 2 F_T(n; k+1, \theta) - F_T(n-1;
k+1, \theta) - F_T(n+1; k+1,\theta) \Big) \qquad n = 0, 1, \dots.
\end{aligned} \\

Which was also found by Cori *et al*^(\[[2](#ref-cori2013new)\]).

### 3.4 Log-Normal distribution

The Log-Normal distribution has the density function:

\\ \begin{aligned} f_T(z;\mu, \sigma) &= {1 \over z \sigma \sqrt{2\pi}}
\exp\left( - {(\log(z) - \mu)^2 \over 2 \sigma^2} \right)\\ F_T(z;\mu,
\sigma) &= 0, \qquad z \< 0. \end{aligned} \\ And distribution function:

\\ F_T(z;\mu, \sigma) = \Phi\left( {\log(z) - \mu \over \sigma} \right).
\\ Where \\\Phi\\ is the standard normal distribution function.

#### 3.4.1 Log-Normal partial expectation

We know that the full expectation of the Log-Normal distribution is
\\\mathbb{E}\[T\] = e^{\mu + \frac{1}{2} \sigma^2}\\, which can be
calculated by integration with the integration substitution \\y = (\ln
z - \mu) / \sigma\\. This has transformation Jacobian:

\\ \frac{dz}{dy} = \sigma e^{\sigma y + \mu}. \\

Doing the same integral for the partial expectation, and using the same
integration substitution gives:

\\ \begin{aligned} \int_t^{t+w_P} z~ f_T(z; \mu, \sigma) dz &= {1 \over
\sigma \sqrt{2\pi}} \int_t^{t+w_P} \mathbb{1}(z \geq 0) \exp\left( -
{(\log(z) - \mu)^2 \over 2 \sigma^2} \right) dz \\ &= {1 \over
\sqrt{2\pi}} \int\_{(\ln t - \mu)/\sigma}^{(\ln(t+w_P) - \mu)/\sigma}
e^{\sigma y + \mu} e^{-y^2/2} dy\\ &= {e^{\mu + \frac{1}{2} \sigma^2}
\over \sqrt{2 \pi} } \int\_{(\ln t - \mu)/\sigma}^{(\ln(t+w_P) -
\mu)/\sigma} e^{-(y- \sigma)^2/2} dy \\ &= e^{\mu + \frac{1}{2}
\sigma^2} \Big\[\Phi\Big({\ln(t+w_P) - \mu \over \sigma} - \sigma\Big) -
\Phi\Big({\ln(t) - \mu \over \sigma} - \sigma\Big) \Big\]\\ &= e^{\mu +
\frac{1}{2} \sigma^2} \Delta\_{w_P}F_T(t; \mu + \sigma^2, \sigma).
\end{aligned} \tag{3.6} \\

#### 3.4.2 Survival function of \\S\_{+}\\ for Log-Normal distribution

By substituting equation [(3.6)](#eq:lognormpartexp) into equation
[(3.3)](#eq:disccensunifprim) we can solve for the survival function of
\\S\_+\\ in terms of analytically available functions:

\\ Q\_{S+}(t ;\mu, \sigma) = Q_T(t + w_P;\mu, \sigma) + { 1 \over w_P}
\Big\[ e^{\mu + \frac{1}{2} \sigma^2} \Delta\_{w_P}F_T(t; \mu +
\sigma^2, \sigma) - t\Delta\_{w_P}F_T(t; \mu, \sigma) \Big\] \\

#### 3.4.3 Log-Normal discrete censored delay distribution

By substituting [(3.6)](#eq:lognormpartexp) into
[(3.3)](#eq:disccensunifprim) we get the discrete censored delay
distribution in terms of analytically available functions:

\\ \begin{aligned} f_n &= (n+1) F_T(n+1; \mu, \sigma) + (n-1) F_T(n-1;
\mu, \sigma) - 2 n F_T(n; \mu, \sigma) \\ &- e^{\mu + \frac{1}{2}
\sigma^2} \Delta_1^{(2)}F_T(n-1;\mu + \sigma^2, \sigma) \\ &= (n+1)
F_T(n+1; \mu, \sigma) + (n-1) F_T(n-1; \mu, \sigma) - 2 n F_T(n; \mu,
\sigma) \\ &+ e^{\mu + \frac{1}{2} \sigma^2} \Big( 2 F_T(n; \mu +
\sigma^2, \sigma) - F_T(n+1; \mu + \sigma^2, \sigma) - F_T(n-1; \mu +
\sigma^2, \sigma) \Big)\qquad n = 0, 1, \dots. \end{aligned} \\

### 3.5 Weibull distribution

The Weibull distribution has the density function:

\\ f_T(z;\lambda,k) = \begin{cases}
\frac{k}{\lambda}\left(\frac{z}{\lambda}\right)^{k-1}e^{-(z/\lambda)^{k}},
& z\geq0 ,\\ 0, & z\<0, \end{cases} \\ And distribution function:

\\ F_T(z;\lambda,k))=\begin{cases}1 - e^{-(z/\lambda)^k}, & z\geq0,\\ 0,
& z\<0.\end{cases} \\ Where \\\Phi\\ is the standard normal distribution
function.

#### 3.5.1 Weibull partial expectation

We know that the full expectation of the Weibull distribution is
\\\mathbb{E}\[T\] = \lambda \Gamma(1 + 1/k)\\, which can be calculated
by integration using the integration substitution \\y = (z /
\lambda)^k\\, which has transformation Jacobian:

\\ \frac{dz}{dy} = \frac{\lambda}{k}y^{1/k - 1}. \\

Doing the same integral for the partial expectation, and using the same
integration substitution gives:

\\ \begin{aligned} \int\_{t}^{t+w_P} z~ f_T(z; \lambda,k) dz &=
\int_t^{t+w_P} \mathbb{1}(z \geq 0)
\frac{kz}{\lambda}\left(\frac{z}{\lambda}\right)^{k-1}e^{-(z/\lambda)^{k}}
dz \\ &= k\int_t^{t+w_P} \mathbb{1}(z \geq 0)
\left(\frac{z}{\lambda}\right)^{k}e^{-(z/\lambda)^{k}} dz \\ &= \lambda
k \int\_{(t / \lambda)^k}^{((t + w_P) / \lambda)^k} \mathbb{1}(y \geq 0)
y y^{1/k - 1} e^{-y} dy \\ &= \lambda\int\_{(t / \lambda)^k}^{((t + w_P)
/ \lambda)^k} \mathbb{1}(y \geq 0) y^{1/k} e^{-y} dy\\ &= \lambda
\Delta\_{w_P} g(t; \lambda,k) \end{aligned} \tag{3.7} \\

Where

\\ g(t; \lambda, k) = \gamma\left(1 + 1/k, \left({t\vee 0 \over
\lambda}\right)^k\right) = \frac{1}{k}\gamma\left(1/k, \left({t\vee 0
\over \lambda}\right)^k\right) -
\frac{t}{\lambda}\exp\left(-\left({t\vee 0 \over
\lambda}\right)^k\right) \\ is a reparametrisation of the lower
incomplete gamma function. Note that the \\\vee\\ operator \\t \vee 0 =
\text{max}(0, t)\\ comes into the expression due to \\\mathbb{1}(y \geq
0)\\ term in the integrand.

### 3.6 Survival function of \\S\_{+}\\ for Weibull distribution

By substituting equation [(3.7)](#eq:weibullpartexp) into equation
[(3.3)](#eq:disccensunifprim) we can solve for the survival function of
\\S\_+\\ in terms of analytically available functions:

\\ Q\_{S+}(t ;\lambda,k) = Q_T(t + w_P; \lambda,k) + { 1 \over w_P}
\Big\[ \lambda \Delta\_{w_P} g(t; \lambda,k) - t\Delta\_{w_P}F_T(t;
\lambda,k)\Big\]. \\

#### 3.6.1 Weibull discrete censored delay distribution

By substituting [(3.7)](#eq:weibullpartexp) into
[(3.3)](#eq:disccensunifprim) we get the discrete censored delay
distribution in terms of analytically available functions:

\\ \begin{aligned} f_n &= (n+1)F_T(n+1) + (n-1)F_T(n-1) - 2nF_T(n) -
\Delta_1\Big\[ \int\_{n-1}^n f_T(z) z ~dz \Big\] \\ &= (n+1)F_T(n+1) +
(n-1)F_T(n-1) - 2nF_T(n) - \lambda \Delta_1^{(2)} g(n-1; \lambda,k) \\
&= (n+1)F_T(n+1) + (n-1)F_T(n-1) - 2nF_T(n) \\ &+ \lambda \[2 g(n;
\lambda,k) - g(n+1; \lambda,k) - g(n-1; \lambda,k)\] \qquad n = 0, 1,
\dots. \end{aligned} \\

Which was also found by Cori *et al*^(\[[2](#ref-cori2013new)\]).

## 4 Learning more

- For more mathematical background on the analytic solutions see the
  [`vignette("why-it-works")`](https://primarycensored.epinowcast.org/dev/articles/why-it-works.md).
- For a more introductory explanation of the primary event censored
  distribution see the
  [`vignette("primarycensored")`](https://primarycensored.epinowcast.org/dev/articles/primarycensored.md).

## References

1\.

Park, S. W., Akhmetzhanov, A. R., Charniga, K., Cori, A., Davies, N. G.,
Dushoff, J., Funk, S., Gostic, K., Grenfell, B., Linton, N. M.,
Lipsitch, M., Lison, A., Overton, C. E., Ward, T., & Abbott, S. (2024).
Estimating epidemiological delay distributions for infectious diseases.
*bioRxiv*. <https://doi.org/10.1101/2024.01.12.24301247>

2\.

Cori, A., Ferguson, N. M., Fraser, C., & Cauchemez, S. (2013). A new
framework and software to estimate time-varying reproduction numbers
during epidemics. *American Journal of Epidemiology*, *178*(9),
1505–1512.
