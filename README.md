## IMR-based sequential BOED

![overview](overview.png)

Given a modeling parameter, $`\mathbf{\theta}=\{\mathcal{M},\, \mathbf{\theta}_{\mathcal{M}}\}`$, which includes a constitutive model and its material properties, and a design $\mathbf{d}$ that describes the experimental setup (e.g., the equilibrium radius), the Inertial Microcavitation Rheometry (IMR) approach numerically solved the spherically symmetric motion of bubble dynamics. 
* The current numerical implementation is based on IMR V1. For more details, see [here](https://github.com/InertialMicrocavitationRheometry/IMR_v1).

In computation, the complete flow states $\mathbf{q}$ include bubble radius, bubble-wall velocity, temperature, and other variables, but they are only partially observable and are denoted as $\mathbf{y}$. The process for the IMR-based sequential BOED is described as follows:
* Optimal design: Maximizing the expected information gain (EIG) using Bayesian optimization (BO) to design the most informative cavitation experiments.
* Model inference:
  * Data assimilation: Characterizing the unknown material properties, $\mathbf{\theta}_{\mathcal{M}}$, by analyzing the bubble dynamics trajectories, $\mathbf{y}$, using En4D-Var.
  * Bayesian model selection: Using the marginal likelihood to calibrate the probability of each constitutive model, $\mathcal{M}$.
    
When the prior is updated using the posterior, one iteration of the sequential design is completed. Soft material properties are shown to be accurately and efficiently characterized by iterating optimal design and model inference processes.
