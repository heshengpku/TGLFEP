# TGLFEP
`TGLF` code called for EP tansport simulations

---

Three input-parameter modules (`mode_in` in `TGLFEP_tglf_map.f90`) are 
`GA-std case`, 
`DIII-D NBI case` (with GYRO inputs), 
`EPtran inputs`.

Different working modes (`code_in` in `TGLFEP_mainsub.f90`): 

1. `TGLFEP_ky` (a single ky) to get growthrate and frequency; (`code_in = 0`)
2. `TGLFEP_TM` (a ky spectrum) to get the growthrate and frequency spectra; 
3. `TGLFEP_ky_widthscan` (width scan at a single ky) to find a suitable `width`; 
4. `TGLFEP_ky_nEPscan` (nEP scan at a single ky) to get the growthrate (and frequency) vs nEP; 
5. Find a suitable `width` and then get the spectrum;
6. `TGLFEP_bisection` to find the density threshold of 
![](http://latex.codecogs.com/gif.latex?\gamma_{AE}>0) with a given `width`;
7. Find a suitable `width` and then find the density threshold of 
![](http://latex.codecogs.com/gif.latex?\gamma_{AE}>0). (`code_in = 6`)

`TGLFEP` can work with `MPI` for
radial scan (`ir`), 
safety factor q scan (`q_factor`), 
shear s scan (`shat_factor`), 
EP density gradient scan (`scan_factor`), 
and the different toroidal n (`n_toroidal`) to find a minimum threshold.

---

> [More introductions on TGLF input](https://fusion.gat.com/theory/Tglfinput)

> The `p_prime` in `TGLF` is usually negative and should use the total pressure including the EP species.
