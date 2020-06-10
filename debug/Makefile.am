AUTOMAKE_OPTIONS = foreign subdir-objects

EMAIL = alexander.sodt@nih.gov

bin_PROGRAMS = betaHist.dbg betaPostProcess.dbg fast_dg.dbg histCollect.dbg genScatteringLengths.dbg bestFitScaling.dbg deuterate.dbg

AM_CPPFLAGS = -include ../config.h $(AM_INCLUDE) 
AM_CFLAGS = -O2 -std=gnu11
AM_CXXFLAGS = -O2

SRC_DIR = ../src

vpath %.C $(SRC_DIR)
vpath %.c $(SRC_DIR)
vpath %.h $(SRC_DIR)

all_SRCS 	= util.C pdb.C alignSet.C input.C hs.C beta_util.C spline.C dcd.C 

betaHist_dbg_SOURCES  		= betaHist.C $(all_SRCS)
betaPostProcess_dbg_SOURCES  	= betaPostProcess.C $(all_SRCS)
deuterate_dbg_SOURCES 		= deuterate.C deuteration.C $(all_SRCS)
fast_dg_dbg_SOURCES   		= fast_dg.c dg_dcd.c bfile.c 
genScatteringLengths_dbg_SOURCES	= genScatteringLengths.C $(all_SRCS)
bestFitScaling_dbg_SOURCES	= bestFitScaling.C NNLS.C $(all_SRCS)
histCollect_dbg_SOURCES		= histCollect.C util.C

betaHist_dbg_LDFLAGS  			= $(AM_LDFLAGS)
betaPostProcess_dbg_LDFLAGS  		= $(AM_LDFLAGS)
deuterate_dbg_LDFLAGS 			= $(AM_LDFLAGS)
fast_dg_dbg_LDFLAGS   			= $(AM_LDFLAGS)
genScatteringLengths_dbg_LDFLAGS	= $(AM_LDFLAGS)
bestFitScaling_dbg_LDFLAGS		= $(AM_LDFLAGS)
histCollect_dbg_LDFLAGS			= $(AM_LDFLAGS)

betaHist_dbg_LDADD  			= $(AM_LDADD) -lstdc++
betaPostProcess_dbg_LDADD  		= $(AM_LDADD) -lstdc++
deuterate_dbg_LDADD 			= $(AM_LDADD) -lstdc++
fast_dg_dbg_LDADD   			= $(AM_LDADD) -lstdc++
genScatteringLengths_dbg_LDADD		= $(AM_LDADD) -lstdc++
bestFitScaling_dbg_LDADD		= $(AM_LDADD) -lstdc++
histCollect_dbg_LDADD			= $(AM_LDADD) -lstdc++



