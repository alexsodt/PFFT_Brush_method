# if your system uses the CHARMM lipid forcefield the programs will likely recognize your atom names.
# for other forcefields you may wish to generate a list of scattering lengths, one per line.

# generates dppc_dfpc.bz [beta(z)], takes a few seconds.
../optimized/betaHist.opt gen_bz.inp dppc_dfpc.psf dppc_dfpc.dcd > gen_bz.out
# runs laterally-averaged
../optimized/betaPostProcess.opt lat.inp > lat.out
# generates distribution of scatterers for the PFFT method. Takes a few minutes.
../optimized/betaHist.opt gen_br.inp dppc_dfpc.psf dppc_dfpc.dcd > gen_br.out
# runs PFFT
../optimized/betaPostProcess.opt pfft.inp > pfft.out
