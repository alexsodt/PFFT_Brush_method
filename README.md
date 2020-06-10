# PFFT_Brush
Software to compute the scattering intensity from a lipid bilayer, including structural correlations in the bilayer plane.

Alpha version written by Mitchell Dorrell and Alex Sodt

Obviously this software is in its early stages and we'd like (even negative) feedback.

There is an example (in the example directory to run).

This will work best for CHARMM forcefields. For PFFT we've made a number of assumptions regarding atom names and scattering lengths.

Two programs are required to generate the intensity.

betaHist.exe: 

In a first pass, here lateral averaging mode, it generates the scattering length density as a function of z (height along bilayer normal).
It takes an input file with a unique name:
e.g.:

name	lo

(For example, see examples/lo/lo_lat_setup.inp)

The syntax is

betaHist.exe input_file psf/pdb dcd

where the second argument can either be a psf (protein structure file) or dcd (NAMD trajectory).

In a second pass, to compute the PFFT intensity, it generates the distribution of scattering pair intensities.

This also requires an input file

name	lo
CUTOFF  30.0

where the value following "CUTOFF" is the PFFT cutoff.
(For example, see examples/lo/lo_pfft_setup.inp)

The second program to use is betaPostProcess.exe, which takes in an input file and produces the intensity.

For example:

hist	lo.br
betaz	lo.bz
name	lo
CUTOFF  30.0

where lo.br is produced by betaHist.exe in PFFT mode, and lo.bz is produced by betaHist.exe in lateral-averaging mode.

email: alexander.sodt@nih.gov
