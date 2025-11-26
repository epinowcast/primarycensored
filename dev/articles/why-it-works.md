# Why it works

## 1 Introduction

### 1.1 What are we going to do in this vignette

In this vignette, we’ll explain the statistical model that the
`primarycensored` package is based on. We’ll cover the following key
points:

1.  Introduction to the censored data problems in time to event
    analysis.
2.  Discuss the relevant issues from censoring in epidemiological data.
3.  Introduce the statistical model used in `primarycensored`.
4.  Distributions where we can derive the censored delay distribution
    analytically.

If you are new to the package, we recommend that you start with the
[`vignette("primarycensored")`](https://primarycensored.epinowcast.org/dev/articles/primarycensored.md)
vignette.

## 2 Censoring and right truncation problems in time to event analysis

Time-to-event analysis, also known as survival analysis, concerns
estimating the distribution of delay times between events. A distinctive
feature of the field are the methodological techniques used to deal with
the missing data problems common in data sets of delay
times^(\[[1](#ref-leung1997censoring)\]). In `primarycensored` we focus
on two particular missing data problems:

- **Interval censoring**. The primary (start) time and/or the secondary
  (end) time of the delay are unobserved but both are known to be within
  intervals.
- **Right truncation**. Truncated delay data items are only reported if
  their secondary events are less than a known time, for example the
  current data collection time.

In statistical epidemiology, these missing data problems occur
frequently in both data analysis and theoretical modelling. For a more
detailed description of these problems in an epidemiological context
see^(\[[2](#ref-Park2024),[3](#ref-Charniga2024)\]).

In data analysis, events in epidemiology are commonly reported as
occurring on a particular day or week (*interval censoring*). In an
emerging outbreak, datasets can be incompletely observed (*right
truncation*) and their can be a great deal of uncertainty around the
precise timing of events (*interval censoring*).

In theoretical epidemiological modelling, it is often appropriate to
model the evolution of an infectious disease as occurring in discrete
time, for example in the [`EpiNow2`](https://epiforecasts.io/EpiNow2/)
and [`EpiEstim`](https://mrc-ide.github.io/EpiEstim/index.html)
modelling packages. This means appropriately discretising continuous
distributions, such as the generation interval distribution. In
`primarycensored` we treat the discretisation of intrinsically
continuous distributions as an *interval censoring* problem which allows
us to simultaneously provide methods for both applied and theoretical
contexts.

## 3 Statistical model used in `primarycensored`

As described in [Getting Started with
`primarycensored`](https://primarycensored.epinowcast.org/dev/articles/primarycensored.md),
`primarycensored` focuses on a subset of methods from time to event
analysis that address data missingness problems commonly found in
epidemiological datasets. We present the statistical problem as a double
interval censoring problem, where both the primary event time and the
secondary event times are interval censored. We can recover single
interval censoring problems by reducing one of the intervals to a point.
In particular, all of the methods in `primarycensored` start by assuming
no secondary event censoring and solving the primary censoring problem.
This is the reason for the packages name despite the most common use
case being when there is double censoring. A key assumption we make
throughout is that the censoring window for events are known and
independent of the event time within the censoring interval. This is
known as non-informative censoring.

The target for inference is the distribution of the delay time between
the primary and secondary events. We assume that the delay time is a
random variable \\T = S - P\_{u}\\ with distribution function \\F_T(t) =
Pr(T \< t)\\ and density function \\f_T(t)\\. In this treatment we
assume that the delay time is shift-invariant, that is, the distribution
of the delay time is the same regardless of the primary event time.

The (unconditional) primary event time is a random variable \\P\_{u}\\
with distribution function \\F\_{P\_{u}}(t)\\ and density function
\\f\_{P\_{u}}(t)\\. The secondary event time is a random variable \\S\\,
but in this treatment we construct the secondary event time from the
primary time and delay, therefore the marginal distribution of \\S\\ is
not considered.

The *censoring window* for each event is the interval within which each
event is known to have occurred, respectively, \\P \in \[t_P, t_P +
w_P\]\\ and \\S \in \[t_S, t_S + w_S\]\\. The lengths of the censoring
windows are \\w_P\\ and \\w_S\\ respectively. The precise event times
within their windows are unobserved. Note that since the primary event
time is known up to the censoring window, we are predominantly
interested in the *conditional* primary time \\P = P\_{u} \| \\ P\_{u}
\in \[t_P, t_P + w_P\]\\\\ which has density function:

\\ \begin{aligned} f\_{P}(p) &= {f\_{P\_{u}}(p) \over F\_{P\_{u}}(t_P +
w_P) - F\_{P\_{u}}(t_P)}, \qquad &p \in \[t_P, t_P + w_P\],\\ f\_{P}(p)
&= 0, \qquad &\text{otherwise}. \end{aligned} \\

In this note, we measure the *censored delay time* \\T_c\\ between the
primary and secondary event windows from endpoint to endpoint: \\T_c =
t_S + w_S - (t_P + w_P)\\. Note that in the [generative model for
delays](https://primarycensored.epinowcast.org/dev/articles/primarycensored.html#generating-random-samples-with-rprimarycensored)
from “Getting started” the truncated delay \\t\_{\text{valid}}\\ is
measured from startpoint to startpoint of event windows.

In our treatment below we focus on the [survival
function](https://en.wikipedia.org/wiki/Survival_function) of the time
after the end of the primary window to the secondary event time which we
denote \\S\_{+}\\. We then use this to derive the distribution of the
censored delay time \\T_c\\. This is equivalent to, but differs in
mathematical approach from other treatments of the censoring problems in
epidemiology, such as^(\[[2](#ref-Park2024)\]), see section [Connections
to other approaches](#connections-to-other-approaches) for details.

### 3.1 Censored delay time distribution

In this section, we explain how to derive the distribution of the
censored delay time \\T_c\\ from the distribution of the delay time
\\T\\ and the condition distribution of the primary event time \\P\\.

#### 3.1.1 Survival function of time from the end of the primary censoring window to the secondary event time

When reasoning about the distribution of the censored delay time
\\T_c\\, it is useful to consider the time from the end (right) point of
the primary censoring interval to the secondary time as a random
variable,

\\ S\_+ = S - (t_P + w_P) = T - ((t_P + w_P) - P) = T - C_P. \\

Where \\T\\ is the delay distribution of interest and \\C_P = (t_P +
w_P) - P\\ is interval between the end (right) point of the primary
censoring window and the primary event time; note that by definition
\\C_P\\ is not observed but we can relate its distribution to the
distribution of \\P\\: \\F\_{C_P}(p) = Pr(C_P \< p) = Pr(P \> t_P +
w_P - p)\\.

With non-informative censoring, it is possible to derive the upper
distribution function of \\S\_+\\, or [*survival
function*](https://en.wikipedia.org/wiki/Survival_function) of \\S\_+\\,
from the distribution of \\T\\ and the distribution of \\C_P\\.

\\ \begin{equation} \begin{split} Q\_{S\_+}(t) &= Pr(S\_+ \> t) \\ &=
Pr(T \> C_P + t) \\ &= \mathbb{E}\_{C_P} \Big\[Q_T(t + C_P)\Big\] \\ &=
\int_0^{w_P} Q_T(t + p) f\_{C_P}(p) dp. \end{split} \end{equation} \\

Using [integration by
parts](https://en.wikipedia.org/wiki/Integration_by_parts) gives: \\
Q\_{S\_+}(t) = Q_T(t + w_P) + \int_0^{w_P} f_T(t+p) F\_{C_P}(p) dp.
\tag{3.1} \\

Where we have used that \\Q^{'}\_{T} = - f_T\\, \\Q_T\\ is the survival
function of the actual delay distribution of interest and \\w_P\\ is the
length of the primary censoring window.

Equation [(3.1)](#eq:survivalfunc) is the key equation in this note and
is used to derive the distribution of the censored delay time \\T_c\\.
It has the interpretation that the probability that the secondary event
time is greater than \\t\\ after the end of the primary censoring window
is the sum of two disjoint event probabilities:

1.  The probability that the *actual* delay time \\T\\ is greater than
    \\t + w_P\\.
2.  The probability that the *actual* delay time \\T\\ is between \\t\\
    and \\t + w_P\\, and the primary event time \\P\\ occurred
    sufficiently close to the end of the primary censoring window that
    the secondary even occurred more than time \\t\\ after the end of
    the primary window.

Note that in [“Getting
started”](https://primarycensored.epinowcast.org/dev/articles/primarycensored.html#compute-the-primary-event-censored-cumulative-distribution-function-cdf-for-delays-with-pprimarycensored)
the target for numerical quadrature is the cumulative distribution
function of the sum of the primary time within the primary censoring
window and the delay time.

#### 3.1.2 Probability of secondary event time within a secondary censoring window

Having constructed the survival function of \\S\_+\\ with equation
[(3.1)](#eq:survivalfunc), using numerical quadrature or in some other
way, we can calculate the probability mass of a secondary event time
falling within a observed secondary censoring window of length \\w_S\\
that begins at time \\n - w_S\\ *after* the primary censoring window.
This is the probability that the censored delay time \\T_c\\ is \\n\\.

This gives the censored delay time probability [by integrating over
censored
values](https://mc-stan.org/docs/stan-users-guide/truncation-censoring.html#integrating-out-censored-values):

\\ Pr(S\_+ \in \[n - w_S, n)) = Q\_{S\_+}(n-w_S) - Q\_{S\_+}(n).
\tag{3.2} \\

Note that the censored secondary event time can also occur within the
primary censoring window. This happens with probability, \\
Q\_{S\_+}(-w_P) - Q\_{S\_+}(0) = 1 - Q_T(w_P) - \int_0^{w_P} f_T(p)
F\_{C_P}(p) dp = Pr(T\< C). \\

## 4 Connections to other approaches

In this section, we discuss how the approach taken in `primarycensored`
relates to some other approaches to the censored data problems in time
to event analysis.

### 4.1 Connection to Park et al 2024

Using the notation from the methods overview given in Park et
al^(\[[2](#ref-Park2024)\]), we write the conditional probability of the
secondary event time \\S\in (S_L,S_R)\\ given the primary event time \\P
\in (P_L,P_R)\\ as:

\\ \begin{aligned} \mathrm{Pr}(S_L \< S \< S_R \| P_L \< P \< P_R) &=
\frac{\mathrm{Pr}(P_L \< P \< P_R, S_L \< S \< S_R)}{\mathrm{Pr}(P_L \<
P \< P_R)} \\ &= \frac{\int\_{P_L}^{P_R} \int\_{S_L}^{S_R} g_P(x)
f_x(y-x) dy dx}{\int\_{P_L}^{P_R} g_P(x) dx}\\ &= \int\_{P_L}^{P_R}
\int\_{S_L}^{S_R} g_P(x\|P_L, P_R) f_x(y-x)dy dx \end{aligned} \\

In this note, we assume that the forward distribution doesn’t vary over
time (such that \\f_x = f\\), then

\\ \int\_{P_L}^{P_R} \int\_{S_L}^{S_R} g_P(x\|P_L, P_R) f_x(y-x)dy dx =
\int\_{P_L}^{P_R} g_P(x\|P_L, P_R) \big\[F(S_R - x) - F(S_L - x)\big\]
dx \\

Then, by using integration by parts, we get:

\\ \begin{split} \int\_{P_L}^{P_R} g_P(x\|P_L, P_R) \big\[F(S_R - x) -
F(S_L - x)\big\] dx &= F(S_R - P_R) - F(S_L - P_R) \\ & -
\int\_{P_L}^{P_R} G_P(x\|P_L, P_R) \big\[f(S_L - x) - f(S_R - x)\big\]
dx \end{split} \tag{4.1} \\ Where we have used that \\\partial_x F(S_R -
x) = - f(S_R - x)\\ and \\\partial_x F(S_L - x) = - f(S_L - x)\\.

We can now compare this to equation [(3.2)](#eq:seccensorprob) by
considering the following transformations:

- \\P_L = -w_P\\ and \\P_R = 0\\, this in this note we are treating the
  endpoint of the primary censoring window as the origin.
- \\S_L = n-w_S\\ and \\S_R = n\\, that is that we are interested in the
  probability of the secondary event time falling within the secondary
  censoring window \\\[n, n+ w_S)\\.

Then equation [(4.1)](#eq:park) becomes:

\\ \begin{aligned} \mathrm{Pr}(S_L \< S \< S_R \| P_L \< P \< P_R) &=
F(n) - F(n-w_S) - \int\_{-w_P}^{0} G_P(x\|-w_P, 0) \big\[f(n - w_S -
x) - f(n - x)\big\] dx \end{aligned} \\ Making the transformation \\x =
-p\\, and rewriting in the notation of this note gives: \\
\begin{aligned} &= F(n) - F(n-w_S) + \int\_{w_P}^{0} G_P(-p\|-w_P, 0)
\big\[f_T(n - w_S + p) - f_T(n +p)\big\] dp \\ &= F(n) - F(n-w_S) +
\int\_{0}^{w_P} G_P(-p\|-w_P, 0) \big\[f_T(n + p) - f_T(n - w_S
+p)\big\] dp\\ &= F(n) - F(n-w_S) + \int\_{0}^{w_P} (1 - F\_{C_P}(p))
\big\[f_T(n + p) - f_T(n - w_S +p)\big\] dp\\ &= F(n + w_P) - F(n-w_S +
w_P) + \int\_{0}^{w_P} \[f_T(n + p - w_S) - f_T(n + p)\] F\_{C_P}(p)
dp\\ &= Q_T(n-w_S + w_P) - Q_T(n + w_P) + \int\_{0}^{w_P} \[f_T(n + p -
w_S) - f_T(n + p)\] F\_{C_P}(p) dp \\ &= Q\_{S\_+}(n-w_S) - Q\_{S\_+}(n
). \end{aligned} \\ which is same as equation
[(3.2)](#eq:seccensorprob).

In this derivation, we have used that \\G_P(x\|-w_P, 0)\\ is the
distribution function from the time *from* the start of the primary
interval *until* primary event time, and \\F\_{C_P}\\ is the
distribution function of the time *until* the end of the primary event
window *from* the primary event time. Therefore, \\G_P(-p\|-w_P, 0) =
Pr(P \< -p \| P \in (-w_P, 0)) = 1 - Pr(C_P \leq p) = 1 - F\_{C_P}(p)\\.

## 5 Learning more

- For more mathematical background on the analytic solutions see the
  [`vignette("analytic-solutions")`](https://primarycensored.epinowcast.org/dev/articles/analytic-solutions.md).
- For a more introductory explanation of the primary event censored
  distribution see the
  [`vignette("primarycensored")`](https://primarycensored.epinowcast.org/dev/articles/primarycensored.md).

## References

1\. Leung, K.-M., Elashoff, R. M., & Afifi, A. A. (1997). Censoring
issues in survival analysis. *Annual Review of Public Health*, *18*(1),
83–104.

2\. Park, S. W., Akhmetzhanov, A. R., Charniga, K., Cori, A., Davies, N.
G., Dushoff, J., Funk, S., Gostic, K., Grenfell, B., Linton, N. M.,
Lipsitch, M., Lison, A., Overton, C. E., Ward, T., & Abbott, S. (2024).
Estimating epidemiological delay distributions for infectious diseases.
*bioRxiv*. <https://doi.org/10.1101/2024.01.12.24301247>

3\. Charniga, K., Park, S. W., Akhmetzhanov, A. R., Cori, A., Dushoff,
J., Funk, S., Gostic, K. M., Linton, N. M., Lison, A., Overton, C. E.,
Pulliam, J. R. C., Ward, T., Cauchemez, S., & Abbott, S. (2024). Best
practices for estimating and reporting epidemiological delay
distributions of infectious diseases. *PLoS Comput. Biol.*, *20*(10),
e1012520. <https://doi.org/10.1371/journal.pcbi.1012520>
