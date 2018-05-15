# This makefile works with the GNU make command, the one find on
# GNU/Linux systems and often called gmake on non-GNU systems, if you
# are using an old style make command, please see the file
# Makefile_oldstyle provided with the package.

# ======================================================================
# Let's start with the declarations
# ======================================================================

# The compiler
FC = gfortran

# flags for debugging or for maximum performance, comment as necessary
#FCFLAGS = -g -fbounds-check
#FCFLAGS = -O2
FCFLAGS = -O2  -Wall -fcheck=all -g -fbacktrace

# the subfolders for modules, source files, object files, data, backup, doc-useful documents
SRC = ./src
MOD = ./mod
OBJ = ./obj

# absolute path for local libraries and modules,  specify for different computer
LIBDIR =  /home/yu/Dropbox/mylib   
LIBMOD =  /home/yu/Dropbox/mylib/mod

# flags to look for all modules (e.g. look for system .mod files, required in gfortran)
FCFLAGS += -I $(MOD)  -I $(LIBMOD) 

# flags to put .mod files for compiled modules 
FMFLAGS = -J $(MOD)  

# libraries needed for linking,  -L declares the absolute path that contains the libraries to be used
# -llib states a library which is lib.a contained in this folder
#LDFLAGS = -li_need_this_lib

 # general libraries, lapack and blas, load by default
LIBGR = -L$(LIBDIR) -llapack -lblas -lgr -lfft
# -lfft 
#FFT =  -L$(LIBDIR)   -lfft
LDEM = -L$(LIBDIR) -lemsys -ljpl -lemat  # used for polar in EM system

# List of executables to be built within the package
#myobj = $(addprefix $(OBJ)/,  $(wildcard *.o) ) 

# update! 20160304 --- For the external subroutines the module calls, they should be complied before the module, or the module will be linked to empty external ones, without complain  and error.
# and use include to split the module into seperate files, avoiding tedious edit work. But it is strongly suggested that save all the related subroutines into one file. 

#dflrtz.o dflrtz_g3.o x2g.o invx2g.o gr_lf.o gr_cjlf.o-- all put into lf_mod 
#inv.o minv.o  ! replaced by LSS solver from Lapack driver
#dfcr_asymt.o  fctn_asymt.o \  ! integrated in gr_po_mod
# 	 dtwofft.o  dfour1.o  mod_twofft.o mod_fftsc.o  fftsc.o \
# 	 gr_rtbp.o  gr_cjrtbp.o \ ! --only used for test

# the second line is to test the 2D torus in PRTBP, put the needed modules at first.
myobj = $(addprefix $(OBJ)/, gr_rk78.o gr_rk78_f90.o  eigrg.o sort_mod.o  eig_sort.o  prt_eigval.o  prt_eigvec.o \
 	  plob.o  deltx.o deltx_lns.o detmat.o  po_mod.o  inv.o  gr_mpower.o  monomat.o  \
 	  poinc_mod.o  poinc.o  poinc_n.o  diffpoinc.o  \
 	  gr_four.o four_seri.o  \
 	  prntft.o  plob_fxd.o fqmax.o  fqext.o  hannin.o \
 	  lf_mod.o  lfunit.o   cj2v_lf.o dvind_dx_lf.o  \
 	  plob_pt4.o eqprop_num.o  fft_mod.o  fft_ob.o  frebas.o  fourier.o   \
 	  rotnum.o  interpol.o  funcs.o  tau_fc.o  gamm.o deriv_gamm.o \
 	  std_map.o \
 	  emconst_mod.o  gr_cjprtbp.o  gr_prtbp.o  deriv_cjprtbp.o  cj2v_prtbp.o \
 	  dvind_dx_prtbp.o  PoincMap_prtbp.o  poinc_tf_prtbp.o  TimeMap_prtbp.o  )
 
lfobj = $(addprefix $(OBJ)/, gr_rk78.o  eigrg.o  sort_mod.o  eig_sort.o  prt_eigval.o  prt_eigvec.o  \
 	   plob.o  plob_y0.o  plob_fxd.o  deltx.o  detmat.o  monomat.o  gr_mpower.o  gr_mfdinst.o \
 	   fft_mod.o fft_ob.o frebas.o  fourier.o  gr_four.o four_seri.o    \
 	   poinc_mod.o poinc.o  poinc_n.o diffpoinc.o     \
 	   ck_mfdinst.o  plob_n.o ck_plob.o  tmap_mfd.o \
 	   plf_mod.o  po_plf_mod.o  PoincMap_tf_plf.o  PoincMap_plf.o  zvc_plf.o \
 	   lf_mod.o  PoincMap_tf_lf.o  PoincMap_lf.o  lfunit.o  \
 	   datan_2pi.o rotnum.o  interpol.o  funcs.o  tau_fc.o  gamm.o gamm_inv.o  deriv_gamm.o \
 	   inv.o null_cal.o  deltx_qr.o curve_time_mod.o  TimeMap_lf.o  dcj_da_debug.o   )
 	   
# 	    to be added later for the correction of the  continuation of the tori
 	 
poobj = $(addprefix $(OBJ)/, gr_rk78.o  plob.o  eigrg.o  sort_mod.o  eig_sort.o   mmeig_sort.o trace.o  prt_eigval.o  prt_eigvec.o \
 	  deltx.o  poinc_mod.o  poinc.o  lfunit.o  inv.o \
 	  lf_mod.o  deltx_lns.o  detmat.o  po_mod.o  monomat.o  \
 	  gr_mpower.o    plob_fxd.o  muller.o)

eigobj = $(addprefix $(OBJ)/,  eigrg.o  sort_mod.o  eig_sort.o mmeig_sort.o trace.o  prt_eigval.o  prt_eigvec.o)



gr_poobj = $(addprefix $(OBJ)/, gr_rk78.o  plob.o  eigrg.o   prt_eigval.o  prt_eigvec.o \
 	     deltx.o  lf_mod.o  cj2v_lf.o  detmat.o  gr_po_mod.o  monomat.o  \
 	     gr_mpower.o    plob_fxd.o)

 
eq_obj = $(addprefix $(OBJ)/,  lf_mod.o ) 
 
testobj = $(addprefix $(OBJ)/, inv.o  deltx.o  null_cal.o )  # for the test of module 
#zvc_plf.o

#myobj = $(addprefix $(OBJ)/, eigrg.o  dflrtz_g3.o prt_eigval.o prt_eigvec.o x2g.o invx2g.o \
#   	 lf_mod.o  po_mod.o  gr_lf.o gr_cjlf.o  dflrtz.o  gr_rtbp.o  gr_cjrtbp.o  \
# 	 monomat.o  gr_mpower.o mfdinst.o lf2cj.o  \
# 	 plob_fxd.o dtwofft.o dfour1.o prntft.o fqmax.o fqext.o)
# 	 
# poinc_n.o poinc_z0.o  --- to be modified
# join into po_mod: pofam.o plob.o fctn.o dfcr.o deltx.o champ.o adams.o  sect.o poinc.o 
#myobj := $(addprefix $(OBJ)/,  $(patsubst %.c,%.o,$(wildcard *.c) ) )

# plot the configuration at different epoches 
mob_lf= $(SRC)/main_ob_lf.f90 
ob_obj =   $(addprefix $(OBJ)/, gr_rk78.o  plob.o  \
 	     lf_mod.o plob_n.o)

ob_lf: $(mob_lf)  $(ob_obj) 
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mob_lf) $(ob_obj) $(LIBGR) -o ob_lf.exe 



# check the unit for lfunit close to the Moon  -- 2018-03-07 10:07:36 
munit= $(SRC)/main_lfunit_moon.f90 
unit_obj = $(addprefix $(OBJ)/,   lfunit_moon.o ) 

unit: $(munit)  $(unit_obj) 
	$(FC) $(FCFLAGS) $(FMFLAGS) $(munit) $(unit_obj) $(LIBGR) -o unit.exe 


# compute the trace of MM 
mtrace= $(SRC)/main_trace.f90
trace: $(mtrace)  $(eigobj) 
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mtrace) $(eigobj) $(LIBGR) -o trace.exe 


# sort pommegv.dat --discard
mmm_eig= $(SRC)/main_mm_eig.f90
mm_eig: $(mmm_eig)  $(eigobj) 
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mmm_eig) $(eigobj) $(LIBGR) -o mm_eig.exe 

# 2D tori in LF 
mlf = $(SRC)/main_curve_lf.f90
curve_lf: $(mlf)  $(lfobj) 
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mlf) $(lfobj) $(LIBGR) -o curve_lf.exe 

mlf_ck = $(SRC)/main_curve_lf_ck.f90
curve_lf_ck: $(mlf)  $(lfobj) 
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mlf_ck) $(lfobj) $(LIBGR) -o curve_lf_ck.exe 

# 2D tori in LF 
mlf_time = $(SRC)/main_curve_time_lf.f90
curve_lf_time: $(mlf_time)  $(lfobj) 
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mlf_time) $(lfobj) $(LIBGR) -o curve_lf_time.exe 


mlf_time_cont = $(SRC)/main_curve_cont_lf.f90
curve_cont_lf: $(mlf_time_cont)  $(lfobj) 
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mlf_time_cont) $(lfobj) $(LIBGR) -o curve_cont_lf.exe 


# 2 curves for check 
mlf2 = $(SRC)/main_2curve_lf.f90
curve_lf2: $(mlf2)  $(lfobj) 
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mlf2) $(lfobj) $(LIBGR) -o curve_lf2.exe 

# starting from a periodic orbit for all the family of tori 
mlf_cont = $(SRC)/main_curve_lf_cont.f90
curve_lf_cont: $(mlf_cont)  $(lfobj) 
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mlf_cont) $(lfobj) $(LIBGR) -o curve_lf_cont.exe 
 
# read from a refined curve for the family  
mcurvehfam_lf = $(SRC)/main_curve_hfam_lf.f90
curve_hfam_lf: $(mcurvehfam_lf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mcurvehfam_lf) $(lfobj) $(LIBGR) -o curve_hfam_lf.exe 
	


# 2D tori in PLF 
mplf = $(SRC)/main_curve_plf.f90
curve_plf: $(mplf)  $(lfobj) 
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mplf) $(lfobj) $(LIBGR) -o curve_plf.exe 

mtori_plf = $(SRC)/main_tori_plf.f90
tori_plf: $(mtori_plf)  $(lfobj) 
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mtori_plf) $(lfobj) $(LIBGR) -o tori_plf.exe 

mtori_fam_plf = $(SRC)/main_tori_fam_plf.f90
tori_fam_plf: $(mtori_fam_plf)  $(lfobj) 
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mtori_fam_plf) $(lfobj) $(LIBGR) -o tori_fam_plf.exe 
	
	
mcurvehfam_plf = $(SRC)/main_curve_hfam_plf.f90
curve_hfam_plf: $(mcurvehfam_plf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mcurvehfam_plf) $(lfobj) $(LIBGR) -o curve_hfam_plf.exe 
	
mcurvefam_plf = $(SRC)/main_curve_fam_plf.f90
curve_fam_plf: $(mcurvefam_plf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mcurvefam_plf) $(lfobj) $(LIBGR) -o curve_fam_plf.exe 



# for the parallel p.o. to x-y plane  
mpo_parall =  $(SRC)/main_po_parall.f90
po_parall:  $(mpo_parall) $(poobj)
	$(FC) $(FCFLAGS) $(mpo_parall) $(poobj)  $(LIBGR) -o po_parall.exe 

#  plot the location of eq 
meq =  $(SRC)/main_eq.f90
eq:  $(meq) $(eq_obj)
	$(FC) $(FCFLAGS) $(meq) $(eq_obj)  $(LIBGR) -o eq.exe 


# asymmetric p.o. 
mpo_asymt =  $(SRC)/main_po_asymt.f90
po_asymt: $(mpo_asymt) $(gr_poobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo_asymt) $(gr_poobj) $(LIBGR) -o po_asymt.exe 


# symmetric periodic computation 
mpo_symt =  $(SRC)/main_po_symt.f90
po_symt: $(mpo_symt) $(poobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo_symt) $(poobj) $(LIBGR) -o po_symt.exe 

mpo_ck =  $(SRC)/main_po_ck.f90
po_ck: $(mpo_ck) $(poobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo_ck) $(poobj) $(LIBGR) -o po_ck.exe 



meq_eig = $(SRC)/main_eqeig.f90
eq_eig: $(meq_eig) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(meq_eig) $(lfobj) $(LIBGR) -o eq_eig.exe 


mmfd_plf = $(SRC)/main_mfd_plf.f90
mfd_plf: $(mmfd_plf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mmfd_plf) $(lfobj) $(LIBGR) -o mfd_plf.exe 


mpomm_plf = $(SRC)/main_pomm_plf.f90
pomm_plf: $(mpomm_plf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpomm_plf) $(lfobj) $(LIBGR) -o pomm_plf.exe 


mzvc_plf = $(SRC)/main_zvc_plf.f90
zvc_plf: $(mzvc_plf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mzvc_plf) $(lfobj) $(LIBGR) -o zvc_plf.exe 


mpcmfd_plf = $(SRC)/main_pcmfd_plf.f90
pcmfd_plf: $(mpcmfd_plf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpcmfd_plf) $(lfobj) $(LIBGR) -o pcmfd_plf.exe 


mtmapmfd_plf = $(SRC)/main_tmapmfd_plf.f90
tmapmfd_plf: $(mtmapmfd_plf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mtmapmfd_plf) $(lfobj) $(LIBGR) -o tmapmfd_plf.exe 


mpomfd_plf = $(SRC)/main_pomfd_plf.f90
pomfd_plf: $(mpomfd_plf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpomfd_plf) $(lfobj) $(LIBGR) -o pomfd_plf.exe 

mpofam_plf = $(SRC)/main_pofam_plf.f90
pofam_plf: $(mpofam_plf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpofam_plf) $(lfobj) $(LIBGR) -o pofam_plf.exe 


mpo_plf = $(SRC)/main_po_plf.f90
po_plf: $(mpo_plf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo_plf) $(lfobj) $(LIBGR) -o po_plf.exe 

mpcob_plf = $(SRC)/main_pcob_plf.f90
pcob_plf: $(mpcob_plf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpcob_plf) $(lfobj) $(LIBGR) -o pcob_plf.exe 


mpoinc_plf = $(SRC)/main_poinc_plf.f90
poinc_plf: $(mpoinc_plf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpoinc_plf) $(lfobj) $(LIBGR) -o poinc_plf.exe 

mfft_lf = $(SRC)/main_fft_lf.f90
fft_lf: $(mfft_lf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mfft_lf) $(lfobj) $(LIBGR) -o fft_lf.exe 


mfft_plf = $(SRC)/main_fft_plf.f90
fft_plf: $(mfft_plf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mfft_plf) $(lfobj) $(LIBGR) -o fft_plf.exe 


mfft_prtbp = $(SRC)/main_fft_prtbp.f90
fft_prtbp: $(mfft_prtbp) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mfft_prtbp) $(myobj) $(LIBGR) -o fft_prtbp.exe 



mmap_plf = $(SRC)/main_map_plf.f90
map_plf: $(mmap_plf) $(lfobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mmap_plf) $(lfobj) $(LIBGR) -o map_plf.exe 

tfur =  $(SRC)/test_fftmod.f 
tfur: $(tfur) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(tfur) $(myobj) $(LIBGR) -o tfur.exe 


mpoinclf = $(SRC)/main_poinc_lf.f90
tori_lf: $(mpoinclf) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpoinclf) $(myobj) $(LIBGR) -o tori_lf.exe 



mtimemap = $(SRC)/main_curve_tmap_prtbp.f90
time_prtbp: $(mtimemap) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mtimemap) $(myobj) $(LIBGR) -o time_prtbp.exe 



mstdmap = $(SRC)/main_curve_stdmap.f90
stdmap: $(mstdmap) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mstdmap) $(myobj) $(LIBGR) -o stdmap.exe 


mtinterpol = $(SRC)/test_interpol.f90
tinterpol: $(mtinterpol) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mtinterpol) $(myobj) $(LIBGR) -o tinterpol.exe 

mcurvehfam_prtbp = $(SRC)/main_curve_hfam_prtbp.f90
curve_hfam_prtbp: $(mcurvehfam_prtbp) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mcurvehfam_prtbp) $(myobj) $(LIBGR) -o curve_hfam_prtbp.exe 
	
mcurvefam_prtbp = $(SRC)/main_curve_fam_prtbp.f90
curve_fam_prtbp: $(mcurvefam_prtbp) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mcurvefam_prtbp) $(myobj) $(LIBGR) -o curve_fam_prtbp.exe 


mprtbp = $(SRC)/main_curve_prtbp.f90
prtbp: $(mprtbp) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mprtbp) $(myobj) $(LIBGR) -o prtbp.exe 


testarea = $(SRC)/test_area.f90
testarea: $(testarea) $(testobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(testarea) $(testobj) $(LIBGR) -o testarea.exe 

mcurvpc =  $(SRC)/main_curve_poinc.f90
curvepc: $(mcurvpc) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mcurvpc) $(myobj) $(LIBGR) -o curvepc.exe 

mrotnum =  $(SRC)/main_rotnum.f90
rotnum: $(mrotnum) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mrotnum) $(myobj) $(LIBGR) -o rotnum.exe 


mpo1 =  $(SRC)/main_po1.f90
mpo2 =  $(SRC)/main_po2.f90

mpo21 =  $(SRC)/main_po21.f90
mpo22 =  $(SRC)/main_po22.f90
mpo23 =  $(SRC)/main_po23.f90
mpo24 =  $(SRC)/main_po24.f90
mpo25 =  $(SRC)/main_po25.f90
mpo26 =  $(SRC)/main_po26.f90
mpo27 =  $(SRC)/main_po27.f90
mpo28 =  $(SRC)/main_po28.f90


mpo3 =  $(SRC)/main_po3.f90
testapo =  $(SRC)/testapo.f90 # use the same data as po3, by call the asymmetric routine to test
testpo_asymt =  $(SRC)/testpo_asymt.f90 # use the same data as po3, by call the asymmetric routine to test

apo22 =  $(SRC)/apo22.f90 # for 2nd family of p.o. with ieq=2

ckhalo =  $(SRC)/ckhalo.f90 # use the same data as po3, by call the asymmetric routine to test


mpots =  $(SRC)/main_pots.f90

mmfd = $(SRC)/main_mfd.f90

meqmf = $(SRC)/main_eqmf.f90

#meq10 =  $(SRC)/main_eq10.f90
mpc10 =  $(SRC)/main_pc10.f90

mpc  =  $(SRC)/main_pc.f90

mfft =  $(SRC)/main_fft.f90
mfftst =  $(SRC)/main_fftst.f90

meq1 =  $(SRC)/main_eq1.f90

mfft1 =  $(SRC)/dxtwofft.f

mtest =  $(SRC)/test.f90

mtestphi =  $(SRC)/testphi.f90

eqformfly =  $(SRC)/eqformfly.f90

poformfly =  $(SRC)/poformfly.f90

mrefncur =  $(SRC)/main_refncur.f90

refncur: $(mrefncur) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mrefncur) $(myobj) $(LIBGR) -o refncur.exe 


mffcur =  $(SRC)/main_fourcur.f90
ffcur: $(mffcur) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mffcur) $(myobj) $(LIBGR) -o ffcur.exe 


mcurve =  $(SRC)/main_curv.f90
curve: $(mcurve) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mcurve) $(myobj) $(LIBGR) -o curve.exe 

mtori =  $(SRC)/main_tori.f90
tori: $(mtori) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mtori) $(myobj) $(LIBGR) -o tori.exe 


poformfly: $(poformly) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(poformfly) $(myobj) $(LIBGR) -o poformfly.exe 



eqformfly: $(eqformly) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(eqformfly) $(myobj) $(LIBGR) -o eqformfly.exe 


testphi: $(mtestphi) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mtestphi) $(myobj) $(LIBGR) -o testphi.exe 


test: $(mtest) $(testobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mtest) $(testobj)  $(LIBGR) -o test.exe 


pots: $(mpots) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpots) $(myobj) $(LIBGR) -o pots.exe 

po1: $(mpo1) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo1) $(myobj) $(LIBGR) -o po1.exe 
	
po2: $(mpo2) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo2) $(myobj) $(LIBGR) -o po2.exe 
	
	
	
# to test the influence of the direction of the initial velocity, we test all the 8 cases, with \pm vz as one group, shown in one window
po21: $(mpo21) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo21) $(myobj) $(LIBGR) -o po21.exe	
po22: $(mpo22) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo22) $(myobj) $(LIBGR) -o po22.exe	
po23: $(mpo23) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo23) $(myobj) $(LIBGR) -o po23.exe	
po24: $(mpo24) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo24) $(myobj) $(LIBGR) -o po24.exe	
po25: $(mpo25) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo25) $(myobj) $(LIBGR) -o po25.exe	
po26: $(mpo26) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo26) $(myobj) $(LIBGR) -o po26.exe	
po27: $(mpo27) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo27) $(myobj) $(LIBGR) -o po27.exe	
po28: $(mpo28) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo28) $(myobj) $(LIBGR) -o po28.exe	
 

po3: $(mpo3) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpo3) $(myobj) $(LIBGR) -o po3.exe 
testapo: $(testapo) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(testapo) $(myobj) $(LIBGR) -o testapo.exe 

testpo_asymt: $(testpo_asymt) $(myobj)  
	$(FC) $(FCFLAGS) $(FMFLAGS) $(testpo_asymt) $(myobj) $(LIBGR) -o testpo_asymt.exe 

apo22: $(apo22) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(apo22) $(myobj) $(LIBGR) -o apo22.exe 


ckhalo: $(ckhalo) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(ckhalo) $(myobj) $(LIBGR) -o ckhalo.exe 



mfd: $(mmfd) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mmfd) $(myobj) $(LIBGR) -o mfd.exe 

# eqmf  eqmf3 for bt=1 and bt=2 
eqmf: $(meqmf) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(meqmf) $(myobj) $(LIBGR) -o eqmf.exe 

pc: $(mpc) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpc) $(myobj) $(LIBGR) -o pc.exe 
 
pc10: $(mpc10) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mpc10) $(myobj) $(LIBGR) -o pc10.exe 

fft: $(mfft) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mfft) $(myobj) $(LIBGR)  -o fft.exe 


fftst: $(mfftst) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mfftst) $(myobj) $(LIBGR) -o fftst.exe 

eq1: $(meq1) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(meq1) $(myobj) $(LIBGR) -o eq1.exe 


fft1: $(mfft1) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(mfft1) $(myobj) $(LIBGR) -o fft1.exe 



pofft =  $(SRC)/fftpo1.f90
pofft: $(pofft) $(myobj)
	$(FC) $(FCFLAGS) $(FMFLAGS) $(pofft) $(myobj) $(LIBGR) -o pofft.exe 



# "make" builds all
#all: $(PROGRAMS)



   
# Another option is to put makefile in the SRC folder, but it then become complicated
# to transfer between different directories, for make and execute command
 
          
   
#apo: $(SRC)/main_apo.f90  $(myprogs) $(basicobj) 
#	$(FC) $(FCFLAGS) $(FMFLAGS) $(SRC)/main_apo.f90 $(myprogs) $(basicobj)  -o apo.exe 



#$intro
# ======================================================================
# Here comes the most interesting part: the rules for prog1, prog2,
# prog3 and prog4, modify to suit your needs
# ======================================================================

# In order to understand the next section, the process of building an
# executable has to be clear: this is typically done in two steps:

# 1.Compilation: every source file required for our program (tipically
# x.f90 or x.F90 in case of Fortran) is compiled into an object file
# (usually x.o)
# 2.Linking: the final executable file is built by "linking" together
# all the object files compiled in the previous step; in this step
# additional pre-compiled libraries of can be added, but this will not
# be treated here.

# These two steps are often performed through the same command and can
# be combined into a single operation, so it is easy to confuse them,
# but, in order to understand what comes further, one has to keep in
# mind that they are logically different operations.


#dependence
#prog2.o: prog2.incf
#prog3: aux.o
#prog4.o: mod.o
#prog4: mod.o

#$conclusion
# ======================================================================
# And now the general rules, these should not require modification
# ======================================================================

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(FC) $(FCFLAGS) $(FMFLAGS) -o $@ -c $^ $(LAPACK)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
$(OBJ)/%.o: $(SRC)/%.f90
	$(FC) $(FCFLAGS) $(FMFLAGS) -o $@ -c $< 

$(OBJ)/%.o: $(SRC)/%.f
	$(FC) $(FCFLAGS) $(FMFLAGS) -o $@ -c $<


$(OBJ)/%.o: $(SRC)/%.F90
	$(FC) $(FCFLAGS) $(FMFLAGS) -o $@ -c $<
	

#$(OBJ)/%.o: $(SRC)/%.f90
#	$(FC) $(FCFLAGS) -c $<

#$(OBJ)/%.o: $(SRC)/%.f
#	$(FC) $(FCFLAGS) -c $<


#$(OBJ)/%.o: $(SRC)/%.F90
#	$(FC) $(FCFLAGS) -c $<
	

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f $(OBJ)/*.o  $(MOD)/*.mod  *.MOD *.o *.mod
first:
	mkdir src obj mod fig dat doc bak


veryclean: clean
	rm -f *~ $(PROGRAMS)


