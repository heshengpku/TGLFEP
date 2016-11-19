# TGLFEP
`TGLF` code called for EP simulations

---

`TGLFEP` code is a program 
to run `TGLF` with EP conveniently, 
where uses `TGLF` (in `GAcode`) 
with the version number 
(on May 7th, 2016) 
''TGLF stable\_r6.0.0-31-gba49'' 
on NERSC ''EDISON\_CRAY Linux x86\_64'' 
or '' CORI Linux x86\_64''

`TGLFEP` source files list as below
* TGLFEP_interface.f90
* TGLFEP_tglf_map.f90
* TGLFEP_ky.f90
* TGLFEP_TM.f90
* TGLFEP_ky_widthscan.f90
* TGLFEP_ky_nEPscan.f90
* TGLFEP_scalefactor.f90
* TGLFEP_mainsub.f90
* TGLFEP_driver.f90

`TGLFEP_interface.f90` contains 
the module of internal parameters
for `TGLFEP` code.

`TGLFEP_tglf_map.f90` contains 
the module of input profiles 
and the subroutine to
read plasma and geometry profiles
from input file `input.profile`
and then map them to `TGLF` 
as the required format.

## Input
Four input profiles (`input.profile`) have been tested: 
`GA-std case`, 
`DIII-D NBI case` (with GYRO inputs), 
`EPtran profiles`,
`Two EP species case`.

## Process
The input parameter `process_in` is
used in `mainsub` subroutine 
to clarify several different simulation scheme:
* `0:` `TGLFEP_ky` to get the modes' wavefunction 
at a given ky (`input.ky`);
* `1:` `TGLFEP_TM` to get growth rate and frequency spectra 
with a ky spectrum 
(![](http://latex.codecogs.com/gif.latex? 0 < k_y \\leq 0.15,\\Delta k_y = 0.01));
* `2:` `TGLFEP_ky_nEPscan` to get growth rate and frequency 
versus ![](http://latex.codecogs.com/gif.latex? n_{EP}) 
at a single ky;
* `3:` (**Recommend**) to get growth rate and frequency spectra 
(and/or the growth rate versus ![](http://latex.codecogs.com/gif.latex? n_{EP})). 
It can run with a given `WIDTH` 
or find a suitable `WIDTH` in a given range 
(set in `input.TGLFEP`);
* `4:` (**Recommend**) to find the density threshold of 
![](http://latex.codecogs.com/gif.latex? \\gamma_{AE}=0) 
for multiple toroidal n (equivalently ky) 
to get a minimum density threshold 
or for a series ![](http://latex.codecogs.com/gif.latex? a/L_{n_{EP}}) 
to calculate ![](http://latex.codecogs.com/gif.latex? C_R), 
with a given `WIDTH` 
or find a suitable `WIDTH` in a given range 
(set in `input.TGLFEP`).

## Parametric scan
`TGLFEP` use `MPI` to do
radial scan (`ir`), 
safety factor q scan (`q_factor`), 
shear s scan (`shat_factor`), 
EP density gradient scan (`scan_factor`), 
or other parametric scan 
but needed to revise the source file `TGLFEP_tglf_map.f90`.

---

> [More introductions on TGLF input](https://fusion.gat.com/theory/Tglfinput)

> The `p_prime` in `TGLF` is usually negative and should use the total pressure including the EP species.
