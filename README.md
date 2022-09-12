# Ising Model implementation in C++

The above code uses the Metropolis–Hastings MCMC algorithm  to simulate the 2D Ising Model. We first perform a series of measurements of various observables, namely - $E, E^2, |M|, M^2, M^4, M$, where $E$ is the energy density and $M$ is the magnetization density in the system. On collecting these measurements, we calculate the autocorrelation times of all the observables (except $M$), and define the maximum of them as the autocorrelation time.

Once we obtain the autocorrelation time $\tau$, we run the Metropolis algorithm again for $2\tau \times N_\text{out}$ times where the $N_\text{out}$ is our desired number of measurements. We now take measurements of the observables every $2\tau$ MC steps to ensure uncorrelated data. 

On obtaining the uncorrelated measurements for the observales, we calculate the errors in the observables and the derived physical quantities ($\chi, C_V, U_L$) using the Jackknife method. 

As of now, we have used a very naive version of the data-collapse method to evaluate the **critical exponents** without taking into account the Jackknife errors evaluated for our physical quantities. Anders Sandvik's document on [Computational Studies of Quantum Spin Systems](https://arxiv.org/abs/1101.3281) is an excellent resource on all things related to spin system simulations, including finite size scaling analysis and error analysis as well.
