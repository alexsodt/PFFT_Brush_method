AUTOMAKE_OPTIONS = foreign subdir-objects

EMAIL = alexander.sodt@nih.gov

bin_PROGRAMS = betaHist.opt betaPostProcess.opt fast_dg.opt histCollect.opt genScatteringLengths.opt bestFitScaling.opt deuterate.opt

AM_CPPFLAGS = -include ../config.h $(AM_INCLUDE) 
AM_CFLAGS = -O2
AM_CXXFLAGS = -O2

SRC_DIR = ../src

vpath %.C $(SRC_DIR)
vpath %.c $(SRC_DIR)
vpath %.h $(SRC_DIR)

all_SRCS 	= util.C pdb.C alignSet.C input.C hs.C beta_util.C spline.C dcd.C 

betaHist_opt_SOURCES  		= betaHist.C $(all_SRCS)
betaPostProcess_opt_SOURCES  	= betaPostProcess.C $(all_SRCS)
deuterate_opt_SOURCES 		= deuterate.C deuteration.C $(all_SRCS)
fast_dg_opt_SOURCES   		= fast_dg.c dg_dcd.c bfile.c 
genScatteringLengths_opt_SOURCES	= genScatteringLengths.C $(all_SRCS)
bestFitScaling_opt_SOURCES	= bestFitScaling.C NNLS.C $(all_SRCS)
histCollect_opt_SOURCES		= histCollect.C util.C

betaHist_opt_LDFLAGS  			= $(AM_LDFLAGS)
betaPostProcess_opt_LDFLAGS  		= $(AM_LDFLAGS)
deuterate_opt_LDFLAGS 			= $(AM_LDFLAGS)
fast_dg_opt_LDFLAGS   			= $(AM_LDFLAGS)
genScatteringLengths_opt_LDFLAGS	= $(AM_LDFLAGS)
bestFitScaling_opt_LDFLAGS		= $(AM_LDFLAGS)
histCollect_opt_LDFLAGS			= $(AM_LDFLAGS)

betaHist_opt_LDADD  			= $(AM_LDADD) -lstdc++
betaPostProcess_opt_LDADD  		= $(AM_LDADD) -lstdc++
deuterate_opt_LDADD 			= $(AM_LDADD) -lstdc++
fast_dg_opt_LDADD   			= $(AM_LDADD) -lstdc++
genScatteringLengths_opt_LDADD		= $(AM_LDADD) -lstdc++
bestFitScaling_opt_LDADD		= $(AM_LDADD) -lstdc++
histCollect_opt_LDADD			= $(AM_LDADD) -lstdc++



