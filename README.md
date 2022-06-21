# Ising Model implementation in C++

To run the plotting scripts for data that is saved in outputs/, following are the instructions - 
1. For ``out_plot.gpl`` file, set the value of $i$ to 4, 6, 12, 20, 22, 25, or 26 for <E/N>, <m>, <|m|>, χ, C_v, U, and τ respectively.
e.g. - ``i=26; load "scripts/out_plot.gpl"``
2. For ``cor_out.gpl`` file, set the value of $j$ to anything from 0 to 50. The value of $j$ is related to temperature by $T_j = 2 + j\cdot 0.01$. This will generate all autocorrelation plots till that temperature.
e.g. - ``j=15; load "scripts/out_plot.gpl"``
3. For ``measurement_plot.gpl`` file, set the value of $i$ to to anything from 0 to 50. The value of i is related to temperature by $T = 2 + i\cdot 0.01$. This will plot the measurements of |m| at $T_i$ when we are calculating autocorrelations.
e.g. - ``i=10; load "scripts/measurement_plot.gpl"``
