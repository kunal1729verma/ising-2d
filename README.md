# Ising Model implementation in C++

The above code uses the Metropolis–Hastings MCMC algorithm  to simulate the 2D Ising Model. We first perform a series of measurements of various observables, namely - $E, E^2, |M|, M^2, M^4, M$, where $E$ is the energy density and $M$ is the magnetization density in the system. On collecting these measurements, we calculate the autocorrelation times of all the observables (except $M$), and define the maximum of them as the autocorrelation time $\tau$. 

Once we obtain the autocorrelation time $\tau$, we run the Metropolis algorithm again for $2\tau \times N_\text{out}$ times where the $N_\text{out}$ is our desired number of measurements. We now take measurements of the observables every $2\tau$ MC steps to ensure uncorrelated data. 

On obtaining the uncorrelated measurements for the observales, we calculate the errors in the observables and the derived physical quantities ($\chi, C_V, U_L$) using the Jackknife method. 

As of now, we have used a very naive version of the data-collapse method to evaluate the **critical exponents** without taking into account the Jackknife errors evaluated for our physical quantities. We have used the data of $L = 8, 16, 32$ to do the same. A method to tackle this is proposed in https://www.iopb.res.in/~somen/SMB_papers/smb_seno_jpa34.pdf which I'll be implementing in the next few days. 

Another problem faced by the code currently is the enormously long runtime even for reasonable lattice sizes of $L=64$. I'm also currently looking into fixing this issue as soon as possible.
