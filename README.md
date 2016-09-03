# TGLFEP
TGLF called for EP tansport simulations

Three input-parameter modules are 'GA-std case', 'DIII-D NBI case with GYRO inputs', 'EPtran inputs'.

Different working mode (subroutine) in 'mainsub.f90': 
TGLF_ky (a single ky) to get growthrate and frequency; 
TGLF_TM (a ky spectrum) to get the growthrate and frequency spectra; 
TGLF_ky_widthscan (width scan at a single ky) to find a suitable width; 
TGLF_ky_nEPscan (nEP scan at a single ky) to get the growthrate (and frequency) vs nEP; 
TGLF_bisection to find the threshold of gamma_AE > 0 with a given width
.
It can work with MPI for
radial scan (ir), 
safety factor q scan (q_factor), 
shear s scan (shat_factor), 
EP density gradient scan (scan_factor), 
and the different toroidal n (n_toroidal) to find a minimum threshold.

The p_prime is usually negative and should use the total pressure including the EP species.
