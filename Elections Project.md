# Macroecological Laws in Elections
## Grilli's laws
To give a mathematically rigorous description of Jacopo Grilli's (2020) methodological framework, it is imperative to start from the dynamic core that unifies the three observed laws: the Stochastic Logistic Model (SLM).
### 1. The Dynamical System: Stochastic Logistic Model (SLM)
The central premise is that the relative abundance of a species $i$ in a community, denoted as $x_i$, evolves according to a Langevin stochastic differential equation:
$$\frac{dx_i(t)}{dt} = \frac{x_i(t)}{\tau_i} \left( 1 - \frac{x_i(t)}{K_i} \right) + \sqrt{\frac{\sigma_i}{\tau_i}} x_i(t) \cdot \xi_i(t)$$
Where:
- $\tau_i$: Time scale of return to equilibrium (growth rate).
- $K_i$: Carrying capacity of the species in the environment.
- $\sigma_i$: Intensity of environmental variability (multiplicative noise).
- $\xi_i(t)$: Gaussian white noise with $\langle \xi_i(t) \rangle = 0$ and $\langle \xi_i(t)\xi_i(t') \rangle = \delta(t-t')$.
#### Stationary Distribution
From the associated Fokker-Planck equation, it is derived that the stationary probability distribution for the abundance $x_i$ is a Gamma Distribution:
$$\mathcal{P}(x_i; \alpha_i, \beta_i) = \frac{1}{\Gamma(\alpha_i) \beta_i^{\alpha_i}} x_i^{\alpha_i - 1} e^{-x_i / \beta_i}$$
The distribution parameters are related to the dynamic parameters as follows:
- $\alpha_i = \frac{2}{\sigma_i^2 \tau_i} - 1$ (Shape parameter. Represents the inverse of the volatility).
- $\beta_i = \frac{K_i \sigma_i^2 \tau_i}{2}$ (Scale parameter. Represents the typical size of the fluctuations).
### 2. AFD (Abundance Fluctuation Distribution)
The AFD describes how the abundance of a single entity fluctuates across different spatial samples (urns/municipalities). To demonstrate universality, we must strip the party's intrinsic size from its fluctuation profile. We achieve this by rescaling the variable by its own mean: $z_i = x_i / \langle x_i \rangle$.

Knowing that the theoretical mean of a Gamma distribution is $\langle x_i \rangle = \alpha_i \beta_i$, we apply the theorem of change of variables $P(z_i) = P(x_i) \left| \frac{dx_i}{dz_i} \right|$, where $\frac{dx_i}{dz_i} = \alpha_i \beta_i$. The scale parameter $\beta_i$ completely cancels out, yielding a universal distribution dependent solely on the shape parameter $\alpha_i$:
$$\mathcal{P}(z_i; \alpha_i) = \frac{\alpha_i^{\alpha_i}}{\Gamma(\alpha_i)} z_i^{\alpha_i - 1} e^{-\alpha_i z_i}$$
By definition, the rescaled variable $z_i$ has a mean of $1$ and a variance of $1/\alpha_i$.
We transition to the logarithmic space to study these fluctuations. We define our final empirical variable as:
$$y_i = \ln(z_i) = \ln\left(\frac{x_i}{\langle x_i \rangle}\right)$$
In this logarithmic space, the distribution manifests as an **Exp-Gamma**.
#### Computational Implementation (The SciPy Translation)
To fit the empirical $y_i$ data to the theoretical Exp-Gamma model using computational tools like `scipy.stats.loggamma`, a strict mapping between the physics and the software's parameterization is required. SciPy parameterizes the `loggamma` function with three parameters: `c` (shape), `loc` (location/shift), and `scale`.

To strictly respect the underlying Langevin dynamics, the scale parameter in SciPy must be fixed to 1 (`fscale=1`). If the scale is allowed to vary, the software mathematically alters the underlying variable from $x$ to $x^{1/\text{scale}}$, destroying the physical meaning of the SLM. By fixing the scale to 1, the SciPy output directly translates to our thermodynamic parameters:

1. **Shape (`c`):** Corresponds exactly to our physical $\alpha_i$.
2. **Location (`loc`):** Represents the natural logarithm of the scale of the underlying rescaled Gamma distribution. Since the theoretical scale of $z_i$ is $1/\alpha_i$, the physical coherence of the model can be verified by confirming that empirically, $\text{loc} \approx -\ln(\alpha_i)$.
#### Electoral Data: The Treatment of Zeros
**1. Grilli's Method: The Poisson Sampling Assumption**
In Grilli's (2020) framework, the SLM defines the abundance $x_i$ as a continuous, strictly positive state variable ($x_i > 0$). However, empirical data consists of discrete integer counts of DNA reads ($k_{i,s}$). To bridge this gap, Grilli assumes a Poisson sampling process:
$$k_{i,s} \sim \text{Poisson}(N_s \cdot x_{i,s})$$
Where $N_s$ is the total number of reads in sample $s$ (the sampling depth), and $x_{i,s}$ is the true underlying continuous abundance.

Under this assumption, if the actual abundance $x_{i,s}$ is very low, there is a non-zero statistical probability of observing $k_{i,s} = 0$. In Grilli's ecology, zeros mean a lack of sampling depth. To fit the model without discarding the zeros, Grilli convolves the theoretical Gamma distribution with the Poisson sampling error, resulting in a Gamma-Poisson mixture (Negative Binomial Distribution).

**2. The Electoral Data Problem**
In an electoral system, we count $100\%$ of the valid votes, not a sample. The true strength or probability of a party $x_{i,s}$ is continuous, but voters are discrete integers.

If a minor party has an underlying probability of $0.15\%$ in a specific demographic area, and a polling station (urn) contains only $400$ total voters ($N_s = 400$), the expected number of votes is $\lambda = 0.6$. Because human beings cannot be fractionalized, the actual vote count $k_{i,s}$ collapses to an integer, frequently yielding $0$ votes. This is a **Sampling/Discretization Zero**.

We adopt a strict truncation for two reasons:

- **A. Confounding Structural vs. Stochastic Zeros:** In ecology, it is assumed that a bacterium could theoretically exist anywhere in the lake. In politics, the ecosystem is bounded by rigid regional topologies (e.g., nationalist parties only print ballots in specific provinces). If we feed our raw matrix into a Negative Binomial fit, the model would treat a province where a party did not run as a Poisson sampling failure.
- **B. The Logarithmic Positivity Constraint:** We operate directly in the logarithmic space: $y_i = \ln(x_i/\langle x_i \rangle)$. The logarithm is a strictly continuous multiplicative operator; it cannot process an absolute zero ($-\infty$).

**3. Treatment of Zeros: The Two-Step Filter and the Conditional Mean**
To respect the topological boundaries of the electoral space and the mathematical limits of the logarithmic operator, we implement a strict two-step filtering protocol:
1. **Macro-level Filter (Structural Zeros):** We eliminate any spatial domain (province) where the sum of a party's votes is strictly zero. This removes topological absences from the phase space.
2. **Micro-level Filter (Fluctuation Zeros):** We eliminate specific urns where the party received zero votes _prior_ to calculating the mean abundance.

Consequently, we redefine our theoretical parameter $\langle x_i \rangle$. We do not calculate the absolute mean; we calculate the **conditional mean abundance**, denoted as $\langle x_i \mid x_i > 0 \rangle$.

The Exp-Gamma distribution in our electoral framework are designed to characterize the dynamic demographic fluctuations of political formations exclusively in the spatial domains where they successfully exceed the threshold of integer discretization (1 vote).

The continuous multiplicative nature of the SLM and the logarithmic transformation ($y_i = \ln(z_i)$) cannot process an absolute zero, as it represents $-\infty$ (an absorbing state of thermodynamic extinction). A zero in an urn is often a "sampling zero"—an artifact of finite size discretization where the underlying probability of the party exists, but the total number of voters in that specific urn was too small to materialize it into a whole integer vote.

Therefore, our model estimates the **conditional mean abundance**, denoted mathematically as $\langle x_i \mid x_i > 0 \rangle$. Attempting to circumvent this by adding a microscopic constant (e.g., $x + \epsilon$) is a fatal methodological error. It artificially introduces a massive density peak at $-\infty$ in the logarithmic space, severely corrupting the estimation of the true thermodynamic volatility $\alpha_i$ and rendering the Kolmogorov-Smirnov test meaningless.
#### The Universal AFD and the Null Model
To test the hypothesis that the electoral system behaves identically to an ecological biome, we attempt a "data collapse" to generate a Universal AFD. The premise is that if all parties, regardless of their size, are subjected to the same overarching physical constraint—the geometric limit of the Zero-Sum space (the Simplex) acting as a universal carrying capacity $K$—their rescaled fluctuations $z_{i,s}$ should collapse onto a single Exp-Gamma curve governed by a universal $\alpha$.

Conversely, we define the **Gaussian (Normal) Distribution** as our rigorous Null Model in the logarithmic space. If a party's growth is purely governed by multiplicative noise (Geometric Brownian Motion) without a carrying capacity or restoring force pulling it back to equilibrium, its fluctuations in the linear space will be Lognormal, which translates precisely to a Normal distribution in the logarithmic $y$ space.
![[Universal_AFD.png]]

| Distribution  | Parameters                          | $D_{KS}$ |
| :-----------: | :---------------------------------- | -------- |
|   Exp-Gamma   | $\alpha = 2.4118$, $\beta = 0.4146$ | 0.0378   |
| Normal (Null) | $\mu = -0.2214$, $\sigma = 0.6912$  | 0.0482   |
#### The Breakdown of Universality and the Evolution of Alpha
Upon empirical verification, the hypothesis of strict universality across all elections and parties is rejected.
![[Comparison_Parties_AFD.png]]
The parameters, specifically $\alpha$, fluctuate significantly depending on the historical election year. In fact, mapping the evolution of $\alpha$ over time reveals one of the most profound socio-physical discoveries of the electoral system. From the late 1970s to 2023, the parameter $\alpha$ has shown a consistent upward trend.

![[AFD_alpha_Evolution.png]]

Because the variance of the local fluctuations is $1/\alpha$, an increasing $\alpha$ mathematically dictates a decreasing variance. The system is experiencing a hardening of its ecological niches. In the early democratic period, a party's local volatility was high (low $\alpha$, high spatial variance). In the modern era, the electorate has "crystallized". Driven by intense political polarization, ideological bloc fidelity, and the national homogenization of information via digital media, voters are less susceptible to local environmental noise. The system is structurally more rigid at the microscopic (urn) level than ever before.
#### Goodness of Fit: The Kolmogorov-Smirnov Distance
To evaluate the competition between the SLM (Exp-Gamma) and the Null Model (Normal), we rely on the Kolmogorov-Smirnov statistic ($D_{KS}$), which measures the maximum absolute vertical distance between the empirical and theoretical Cumulative Distribution Functions (CDFs).

For massive sample sizes ($N > 10,000$ polling places), traditional $p$-values become hypersensitive to microscopic deviations, so they are not important. Therefore, the absolute distance $D_{KS}$ is a valid metric. A historical inversion occurs between the years 2000 and 2008, where the Normal distribution outperforms the Exp-Gamma ($D_{Normal} < D_{ExpGamma}$). During these years of extreme bipolarization (hegemony of two major parties), minor parties were crushed out of their ecological equilibrium. Stripped of a stable carrying capacity $K$, their dynamics temporarily reverted to free-floating multiplicative noise, perfectly captured by the Gaussian null model.

![[D_KS_Comparison_AFD.png]]

Election | $\alpha$ | $D_{expg}$ | $D_{norm}$
:-: | :-: | :-: | :-:
1982-10 | 1.476662 | 0.053830 | 0.025968
1986-06 | 1.749540 | 0.032221 | 0.034216
1989-10 | 1.936604 | 0.029416 | 0.035754
1993-06 | 2.054603 | 0.032514 | 0.041854
1996-03 | 2.499118 | 0.035007 | 0.053944
2000-03 | 2.314167 | 0.050785 | 0.046473
2004-03 | 2.798645 | 0.043881 | 0.032386
2008-03 | 2.571624 | 0.051897 | 0.032698
2011-11 | 2.437782 | 0.030429 | 0.053714
2015-12 | 2.996219 | 0.023411 | 0.046147
2016-06 | 3.283535 | 0.025269 | 0.051468
2019-04 | 2.978436 | 0.030147 | 0.065438
2019-11 | 3.102733 | 0.027579 | 0.065183
2023-07 | 3.076761 | 0.037376 | 0.065591
#### Size-Based Segregation: The Dual Ecosystem
The rejection of a purely universal AFD leads to the final topological partition of the ecosystem. The dynamics of a party are fundamentally dictated by its mean mass. When segregating the parties into size classes (let's say: large $>15\%$, medium $5-15\%$, small $<5\%$), two distinct thermodynamic regimes emerge:

1. **Large Parties:** Hegemonic parties strictly follow the **Exp-Gamma** distribution. Because their mean abundance is massive, they frequently collide with the absolute physical limit of the ecosystem (the $100\%$ zero-sum wall). This geometric wall acts as a relentless restoring force, generating the long left tail and the abrupt right wall characteristic of the Exp-Gamma. Their shape is dictated purely by the geometry of the space.
2. **Small Parties:** Minor parties overwhelmingly sometimes adapt better to the **Normal (Gaussian)** distribution. Because they float far below the $100\%$ wall, they lack a structural boundary to deform their symmetry. Furthermore, small parties are highly susceptible to phenomena like "strategic voting" (useful vote). In moments of high polarization, their voters abandon them abruptly for larger blocs. These additive historical shocks break the restorative logic of the SLM, collapsing their sum into Gaussian white noise via the Central Limit Theorem.

**Statistical Note on Small Samples:** When evaluating Size-Based segregation, the number of entities (parties) in a specific class drops drastically (often $N < 15$ parties). For small $N$, the critical KS threshold for significance ($D_{crit} \approx 1.36/\sqrt{N}$) is exceptionally forgiving. A $D_{KS}$ around $0.15$ to $0.20$ is statistically robust, validating the differentiation between the physical regimes of large and small electoral entities.

![[Party_Samples_AFD.png]]

For example, in 2011 Elections:

| $D_{KS}$  |   Small    |   Medium   |   Large    |
| :-------: | :--------: | :--------: | :--------: |
| Exp-Gamma |   0.1016   | **0.0428** | **0.0810** |
|  Normal   | **0.0676** |   0.0475   |   0.1130   |

#### Note: Grilli's Z-Score calculations

Grilli standardizes the logarithmic variable to produce a $z$-score representation: $y_{std} = (y - \langle y \rangle) / \sigma_y$.

While visually convenient for collapsing disparate datasets onto a standardized normal-like scale for comparative plotting, applying this $z$-score transformation to electoral time-series fundamentally destroys the underlying physics. In our framework, the variance of the rescaled fluctuations ($1/\alpha$) is not a nuisance parameter to be normalized away; it is the exact physical signal that measures the volatility and crystallization of the ecosystem. By dividing by the standard deviation $\sigma_y$, one artificially forces the variance of every election to be $1$, erasing the historical evolution of the $\alpha$ parameter and blinding the analysis to the thermodynamic hardening of political niches over time.
![[Grillis_trick.png]]
### 3. MAD (Mean Abundance Distribution)
The MAD examines the macroscopic architecture of the ecosystem by characterizing the variation of mean abundances $\langle x_i \rangle$ among all species (parties) in the community. While the AFD describes how a single entity fluctuates, the MAD describes how the total available carrying capacity is partitioned among the competing entities.
#### Empirical Estimation and the Finite Size Constraint
The extraction of empirical data for the MAD is subjected to severe finite-size limitations. For a given election, the dataset consists simply of the conditional mean vote fraction $\langle x_i \rangle$ for each party $i$.
1. **The Impossibility of Universality:** Unlike the AFD, where fluctuations can be rescaled ($z_i = x_i/\langle x_i \rangle$) to collapse all parties and historical elections into a single universal curve, the MAD represents the static architecture of the system at a specific moment in time. Aggregating all elections into a single MAD is physically invalid because the macro-state of the ecosystem (the number of active parties and the degree of polarization) shifts between electoral cycles. Therefore, the MAD must be fitted and evaluated independently for each election.
2. **The CDF over PDF:** A critical constraint in electoral socio-physics is the low number of macroscopic entities. A single election typically contains only $N \approx 30$ to $50$ consolidated parties. With such a low $N$, constructing a Probability Density Function (PDF) via histograms introduces aggressive binning artifacts. Consequently, the analysis must be performed strictly using the Empirical Cumulative Distribution Function (ECDF), which preserves the exact geometric location of every data point without grouping.
#### The Exact Distribution and Null Models

To understand the generative mechanisms of the electoral system, we test the empirical MAD against three competing macroecological frameworks:

**1. The Exact Distribution: Log-Normal (Multiplicative Growth)**

The paper postulates that the means $\langle x_i \rangle$ across species follow a Log-Normal Distribution:

$$\mathcal{P}(\langle x \rangle; \mu, \sigma_{MAD}) = \frac{1}{\langle x \rangle \sigma_{MAD} \sqrt{2\pi}} \exp \left( -\frac{(\ln \langle x \rangle - \mu)^2}{2\sigma_{MAD}^2} \right)$$

_Physical Meaning:_ The Log-Normal distribution is the mathematical signature of **Gibrat's Law of Proportionate Growth**. It implies that the creation of a political niche is governed by multiplicative network contagion (social media, rallies, word-of-mouth). The growth rate of a party's mean size is proportional to its current size, unconstrained by any central restorative force.

**2. Null Model A: The Gamma Distribution (Neutral Theory)**

In standard ecology (e.g., Hubbell's Neutral Theory of Biodiversity), the MAD is often modeled via a Gamma distribution.

_Physical Meaning:_ While visually similar to the Log-Normal, the Gamma distribution implies the existence of a systemic demographic carrying capacity that pulls the ecosystem toward an equilibrium state. It forces a "thinner" right tail, making the existence of hyper-giant hegemonic parties mathematically less probable than in a Log-Normal regime.

**3. Null Model B: The Beta Distribution (Zero-Sum / Broken Stick)**

Because the electoral space is ultimately bounded by $100\%$ of the votes, one could hypothesize that the mean sizes are distributed according to a Beta distribution, which is strictly bounded between $0$ and $1$.

_Physical Meaning:_ This would imply that the static architecture of the parties is purely a geometric artifact of dividing a finite space (a "broken stick" topology).

![[MAD_Sample_Fits.png]]
#### Validation: The Kolmogorov-Smirnov Distance

To rigorously quantify the goodness-of-fit without the artifacts of data binning, we employ the Kolmogorov-Smirnov (KS) test. When evaluating the $D_{KS}$ values for the MAD, the finite system size ($N \approx 30$ to $50$ parties) dictates the statistical interpretation. The critical threshold to reject a model at a $95\%$ confidence level ($\alpha = 0.05$) is approximated by:

$$D_{crit} \approx \frac{1.36}{\sqrt{N}}$$

For $N = 40$, the critical distance is $D_{crit} \approx 0.215$.

![[MAD_D_KS.png]]
#### Nature of the Electoral System

Empirical evaluation consistently demonstrates that the **Log-Normal distribution** systematically outperforms both the Gamma and Beta distributions across almost all electoral cycles.

This establishes a profound duality in the physics of elections. While the spatial fluctuations of a party at the microscopic level of the urn are rigidly bounded by the Zero-Sum geometry of the ballot box (yielding the Exp-Gamma AFD), the genesis of the party's macroscopic mean size $\langle x_i​\rangle$ operates unconstrained in the national information space. Political ecosystems are **not constrained by local carrying capacities**. Instead, they probably operate as scale-free social networks, where the accumulation of voters is governed by preferential attachment, media amplification, and multiplicative contagion—mechanisms homologous to the proportionate growth dynamics observed in financial markets and wealth distribution.
### 4. Taylor's Law: The Scaling of Fluctuations
Taylor's Law establishes a universal empirical power-law relationship between the spatial variance $\sigma^2(x_i)$ and the spatial mean $\langle x_i \rangle$ of species abundances within an ecosystem:
$$\sigma^2(x_i) = C \cdot \langle x_i \rangle^b$$
Where $C$ is a scaling constant and $b$ is the Taylor exponent. In the logarithmic space, this power law translates into a linear relationship, where the slope $b$ dictates how violently the fluctuations of a political party scale as its voter base grows:
$$\log_{10}(\sigma^2) = \log_{10}(C) + b \log_{10}(\langle x \rangle)$$
#### Standardized Major Axis (SMA) Regression
To empirically calculate the exponent $b$, we must perform a linear regression in the log-log space. However, utilizing a standard Ordinary Least Squares (OLS) regression introduces a critical methodological flaw.

OLS minimizes only the _vertical_ residuals, relying on the strict mathematical assumption that the independent variable (the $X$-axis) is measured with absolute perfection, and all noise resides in the dependent variable (the $Y$-axis). In our electoral ecosystem, this assumption is false. Both the spatial conditional mean ($\langle x_i \rangle$) and the spatial variance ($\sigma^2(x_i)$) are empirical estimators derived from finite samples (urns). Consequently, both axes contain natural sampling noise and discretization error.

When OLS is applied to data where the $X$-axis contains intrinsic error, it suffers from "attenuation bias" (or regression dilution), systematically underestimating the true slope. To correct for this, we employ **Standardized Major Axis (SMA) regression**. By computing the slope as $b_{SMA} = b_{OLS} / R$ (where $R$ is the Pearson correlation coefficient), SMA symmetrically accounts for the stochasticity in both dimensions. This prevents the artificial flattening of the curve and yields the true underlying allometric scaling exponent.
#### Theoretical Expectation from the SLM ($b = 2$)
In Grilli's macroecological framework, the theoretical expectation for the Taylor exponent is strictly $b = 2$, derived directly from the fundamental properties of the Stochastic Logistic Model's stationary Gamma distribution.
For a random variable following a Gamma distribution defined by a shape parameter $\alpha_i$ and a scale parameter $\beta_i$, the moments are exactly:
1. Mean: $\langle x_i \rangle = \alpha_i \beta_i$
2. Variance: $\sigma^2(x_i) = \alpha_i \beta_i^2.$

By isolating the scale parameter as $\beta_i = \langle x_i \rangle / \alpha_i$ and substituting it into the variance equation, we obtain the theoretical relationship:
$$\sigma^2(x_i) = \frac{\langle x_i \rangle^2}{\alpha_i}$$
Assuming that the environmental noise (captured by the inverse volatility parameter $\alpha$) is a systemic constant across all species in a given biome, the variance scales precisely with the square of the mean ($\sigma^2 \propto \langle x \rangle^2$), yielding an exponent of $b = 2$.

A Taylor exponent of exactly $2$ acts as the mathematical signature of **pure multiplicative environmental noise**. It implies that fluctuations are strictly proportional to size; an external systemic shock (such as a financial crisis or a national scandal) affects a _percentage_ of a party's electorate, rather than an absolute, additive number of voters.
#### Empirical Reality: The Suppressed Exponent ($1.6 < b < 1.95$)
Upon empirical evaluation across all elections, the electoral ecosystem reveals a structural deviation from the pure mathematical theory. As observed in the historical data, the empirical Taylor exponent $b$ oscillates between $1.60$ and $1.90$, systematically failing to reach the theoretical $b = 2$ threshold.

To interpret this, we must evaluate the boundary conditions of Taylor's Law:

- **If $b = 1$:** The system is driven by pure additive demographic noise. This is equivalent to a Poisson process, where voters act as entirely independent coin flips without any overarching systemic environmental pressure.
- **If $b = 2$:** The system is driven by pure multiplicative noise, with unconstrained proportional scaling.

Because the electoral exponent resides in the upper-middle of this spectrum ($1.6 < b < 1.95$), we conclude that while the system is heavily dominated by multiplicative forces, it is simultaneously subjected to a **variance suppression mechanism** at macroscopic scales.

An exponent less than 2 indicates that the spatial variance of massive, hegemonic parties grows slower than the SLM predicts. In the physical reality of the ballot box, a minor party might see its vote fraction double from $2\%$ to $4\%$ across different municipalities (exhibiting massive relative variance). However, a giant party with a spatial mean of $40\%$ cannot physically double its vote share to $80\%$ in local samples without colliding violently with the absolute $100\%$ limit of the Zero-Sum geometric space.

Large parties inherently possess both a "structural floor" (a baseline of rigidly loyal partisan voters impervious to environmental noise) and a "demographic ceiling." This confinement actively dampens the extreme spatial fluctuations of massive parties, dragging the global scaling exponent $b$ down from the theoretical $2.0$ to the empirical $\sim 1.7$.
#### Historical Evolution of the Scaling Exponent

![[Taylor_Evolution.png]]

The temporal evolution of the $b$ exponent acts as a macro-seismograph, tracking the structural stability and polarization of the Spanish electoral ecosystem:

1. **1993 - 2008:** During this period, the scaling exponent steadily climbed, peaking at $b \approx 1.88$ in 2008. This era was characterized by the absolute hegemony of two major political blocs. The ecosystem was highly predictable and behaved almost exactly like a pure multiplicative environment. The two systemic giants absorbed almost all environmental fluctuations proportionally, pushing the ecosystem to the edge of the theoretical $b=2$ limit.
2. **2011 - 2015:** The exponent suffers a violent collapse, plunging to $b \approx 1.62$ in 2011 and bottoming out at $b \approx 1.59$ in 2015. This thermodynamic crash coincides precisely with the systemic rupture of the Spanish political landscape (the 15-M movement and the abrupt emergence of new national parties). The sudden influx of mid-sized entities broke the established multiplicative scaling. The spatial variance of the old hegemonic parties collapsed as their traditional voter bases shattered and scattered across the new options, forcing the entire ecosystem into a state characterized by heavy additive noise.
3. **2019 - 2023:** In recent electoral cycles, the exponent shows a slow, secular recovery back toward higher values ($b \approx 1.68$ in 2023). This indicates that the fragmented, multi-party ecosystem has survived its initial chaotic inception and is beginning to crystallize into a new stable multiplicative regime, establishing rigid niches and structural predictability once again.
