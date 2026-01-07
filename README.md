# Analytical and Stochastic Modeling of Radioactive Decay Chains

This repository provides a comprehensive MATLAB framework for simulating radioactive decay kinetics. It implements three distinct computational approaches to solve the coupled differential equations governing nuclear transmutations, supporting complex branching ratios and varying timescales (from microseconds to gigayears).

## 1. Physical & Mathematical Framework

The evolution of a decay chain is governed by the **Law of Radioactive Decay**, where the rate of change of a population $N_i$ depends on its own decay constant $\lambda_i$ and the feeding from its parent(s) $N_{i-1}$.

### A. Recursive Approach (Bateman Equations)
For a linear chain $1 \to 2 \to \dots \to n$, the population of the $n$-th isotope at time $t$ is given by:

$$
N_k(t) = \sum_{n=1}^k C_{kn}e^{-\lambda_nt}, \quad \text{con } 
C_{kn} = 
\begin{cases} 
   C_{k-1,n}\frac{\lambda_{k-1}}{\lambda_k-\lambda_n}, & n \neq k \\
   N_k(0) - \sum_{n=1}^{k-1}C_{kn}, & n=k 
\end{cases}
$$

It can be also applied to chains with branching by considering the branching ratios and correctly identifying the parents. The process is implemented in `DecayChain.m`, this method uses recursive coefficients $C_{i,j}$ to solve for $N(t)$.

### B. Matrix Formalism (Linear Algebra)
To handle complex topologies, we represent the system as:

$$\frac{d\mathbf{N}(t)}{dt} = \mathbf{\Lambda} \mathbf{N}(t)$$

Where $\mathbf{\Lambda}$ is the **Transition Rate Matrix**. The analytical solution is obtained via the **Matrix Exponential**:

$$\mathbf{N}(t) = \exp(\mathbf{\Lambda} t) \mathbf{N}_0$$

Implemented in `DecayMatrix.m` using MATLAB's `expm()`.

### C. Stochastic Simulation (Monte Carlo Method)
Implemented in `DecayMonteCarlo.m`, this solver discretizes time and uses **Binomial** and **Multinomial** distributions to simulate the probabilistic nature of decay. This is particularly useful for observing "quantum noise" in low-activity regimes.

## 2. Detailed Physics Analysis

### A. The $^{99}\text{Mo}/^{99m}\text{Tc}$ Medical Generator
The simulation in `Molibdeno_99.m` models the production of Technetium-99m, the most widely used radioisotope in nuclear medicine. This system is a classic example of **Transient Equilibrium**.

* **Decay Characteristics:** The parent isotope ($^{99}\text{Mo}$) has a half-life of $T_{1/2} \approx 6$ h, while the metastable daughter ($^{99m}\text{Tc}$) has a $T_{1/2} \approx 6$ h.
* **Transient Equilibrium:** Since $\lambda_{parent} < \lambda_{daughter}$, after approximately 4 half-lives of the daughter, the ratio of their activities becomes constant.
* **Branching Efficiency:** The model accounts for the fact that only $87.6%$ of $^{99}\text{Mo}$ decays lead to the metastable $^{99m}\text{Tc}$, while the remaining $12.4%$ decay directly to the ground state $^{99}\text{Tc}$.
* **Clinical Relevance:** The simulation identifies the **Maximum Yield Time** ($t \approx 22.85$ h), which is the optimal moment for "milking" the generator to obtain the highest activity of $^{99m}\text{Tc}$ for diagnostic imaging.

![Mo-99 / Tc-99m Kinetics](Activity_comparison_AnvsMC.png)

### B. Natural Decay Series (Geological Timescales)
The script `Natural_Chains.m` allows for the exploration of the four primordial decay series that shaped the Earth's radiogenic heat and isotopic composition:

| Series | Parent Isotope | Final Stable Nucleus | Main Characteristic |
| :--- | :--- | :--- | :--- |
| **Thorium (4n)** | $^{232}\text{Th}$ | $^{208}\text{Pb}$ | Includes the gaseous $^{220}\text{Rn}$ (Thoron). |
| **Neptunium (4n+1)** | $^{237}\text{Np}$ | $^{209}\text{Bi}$ / $^{205}\text{Tl}$ | Mostly extinct in nature; ends in Bismuth/Thallium. |
| **Uranium (4n+2)** | $^{238}\text{U}$ | $^{206}\text{Pb}$ | Source of indoor $^{222}\text{Rn}$ gas hazards. |
| **Actinium (4n+3)** | $^{235}\text{U}$ | $^{207}\text{Pb}$ | Significant for U-Pb dating in geochronology. |

**Technical Implementation Details:**
* **Numerical Stiffness:** To handle the disparity between $^{238}\text{U}$ ($4.47 \times 10^9$ years) and $^{214}\text{Po}$ ($0.164$ ms), the simulation utilizes a **logarithmic time scale** and for the matrix exponential solvers, the smaller half-lifes are artificially increased to maintain numerical stability.
* **Stochastic Noise Analysis:** The Monte Carlo solver (`DecayMonteCarlo.m`) computes **Residuals** ($\Delta N = N_{mc} - N_{ana}$), visualizing the Poissonian fluctuations that dominate when population sizes are small, a phenomenon crucial for sensitive radiation detection.

## 3. Dependencies
- **MATLAB R2021a** or later.
- **Statistics and Machine Learning Toolbox** (for `mnrnd` and `binornd`).
