# Stochastic-wave-on-sphere
Python code for numerical experiments for a stochastic wave equation on the unit sphere

coeff_rand_field_v6.py: create the initial random field

Outputs:
FracPDE_Nside2048_instance1kap2p3.mat  with kapp1 = 2.3
FracPDE_Nside2048_instance1kap4p1.mat       kapp1 = 4.1
FracPDE_Nside2048_instance1kap4p5.mat       kapp1 = 4.5

ErrorsLtr_wRF_wave.py: compute the average square errors
VarLtr_wave.py: compute the variance 
VarUxy_wave.py: compute the variance Var|U(x,t) - U(y,t)| over 100 samples

waveU0.py: generate the map of the initial random field
waveHomoSoln.py: generate homogeneous soln
waveInhomoSoln.py: generate the inhomogenous solution
                   could be used to create the full solution

waveFullSoln.py: generate the full solution

Outputs:
Initial random fields
FracPDE_Nside2048_instance1kap2p3.png  with kapp1 = 2.3
FracPDE_Nside2048_instance1kap4p1.png       kapp1 = 4.1
FracPDE_Nside2048_instance1kap4p5.png       kapp1 = 4.5

Homogenous solutions
Homo_soln_Lmax600_t_10tau_tau_0p04.png     at tau = 10tau with tau = 0.04, alpha = 0.9, kap1 = 2.3

Inhomogenous solution
  Inhom_soln_Lmax600_t_10tau_tau0p04_alpha1.png  at tau=10tau with tau=0.01 alpha=1.0

plt_soln.py
  plot pre-computed solution

All figures

Initial random fields
FracPDE_Nside2048_instance1kap2p3.png
FracPDE_Nside2048_instance1kap4p1.png
FracPDE_Nside2048_instance1kap4p5.png

Full solutions:
Full_soln_Lmax600_t_10tau_tau0p04_alpha0p9_kap1_2p3_kap2_2p5_a.png
Full_soln_Lmax600_t_10tau_tau0p04_alpha0p9_kap1_2p3_kap2_2p5.png
Full_soln_Lmax600_t_10tau_tau0p04_alpha0p9_kap1_4p1_kap2_2p5.png
Homogeneous soln
Homo_soln_Lmax600_t_10tau_tau_0p04_kap1_4p1.png : kappa1 = 4.1
Homo_soln_Lmax600_t_10tau_tau_0p04.png  : kappa1 = 2.3

Inhomogeneous soln
Inhom_soln_Lmax600_t_10tau_tau0p04_alpha0p9_kap1_2p3_kap2_2p5.png
Inhom_soln_Lmax600_t_10tau_tau0p04_alpha0p9_kap1_4p1_kap2_2p5.png
Inhom_soln_Lmax600_t_10tau_tau0p04_alpha0p9.png
Inhom_soln_Lmax600_t_10tau_tau0p04_alpha1.png
Inhom_soln_Lmax600_t_10tau_tau0p04_instance1.png

Initial condition
Initial_cond_Lmax600_t0_kap1_4p1.png
Initial_cond_Lmax600_t0.png  : kappa1 = 2.3
