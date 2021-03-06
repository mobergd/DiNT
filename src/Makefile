#Makefile for dint

#set POT here or use 'gmake POT=potname'

#POT=dd
#POT=h2oM
#POT=ch4hALL
#POT=ch4ohALL
#POT=ethylene
POT=tbplusexp6-2

EXEDIR = ../exe
TARG = $(EXEDIR)/dint2-$(POT).x$(EXESUFFIX)
OBJDIR = obj
OBJDIR_OPT = obj_opt

FC = mpiifort
LD = mpiifort
OPT_FFLAGS = -O3
FTPPFLAGS = -I../../sprng/SRC -I.. -c
FFLAGS = -c
COMPLIB = -mkl
RNGLIB = -lmlfg
LIBS = $(COMPLIB) $(RNGLIB)
# LDFLAGS = -O3 -g -L../../sprng/lib -L/soft/apps/packages/lapack-3.0/lib
LDFLAGS = -O3 -g -L../../sprng/lib -L/usr/lib64 

VPATH = ..:../../pot:../../dd_pot:../../autopot

SRC = $(POT) \
	addrot.F \
	andersen.F \
	angmom.f \
	ange.f \
	atomdiatom.F \
	bootstrap.f \
	bsstep.f \
	carttojac.f \
	checkfrus.f \
	checkhop.F \
	cn.f \
	decocheck.F \
	derivs.f \
	detwell.f \
	diamin.f \
	diamin2.f \
	diapot.f \
	diapot2.f \
	diatom.f \
	dint.F \
	driver.f \
	ejscale.f \
	elecdeco.f \
	endmpi.f \
	ems.F \
	ewkb.f \
	ewkb2.f \
	finalstate.f \
	getdvec2.f \
	getgrad.f \
	getpem.f \
	getrel.f \
	getrelbz2.f \
	getrho.f \
	gettemp.f \
	header.f \
	hessems.f \
	hop.f \
	honey.f \
	initelec.f \
	initmol.f \
	initmpi.f \
	jointprint.F \
	lindemann.f \
	mmid.f \
	noang.f \
	normod.f \
	nmturn.f \
	nmpot.f \
	orient.f \
	period.f \
	period2.f \
	pjsplit.f \
	getplz.f \
	popnorm.F \
	preatomdiatom.f \
	prediatom.f \
	premol.f \
	preptraj.f \
	pzextr.f \
	radialdist.f \
	ranclu.F \
	rancol.F \
	rangaussian.F \
	ranj.F \
	ranno.F \
	ranrot.F \
	ransamp.F \
	ranstate.F \
	rantherm.F \
	rarray.f \
	readin.f \
	rk4.f \
	rtox.f \
	spin.f \
	stodeco.f \
	takestep.f \
	timing.f \
	torphi.f \
	turn.f \
	turn2.f \
	vibwells.F \
	vwkb.f \
	vwkb2.f \
	xptoy.f \
	airy.f

OBJ = $(patsubst %.f,%.o,$(SRC))
OBJ := $(patsubst %.F,%.o,$(OBJ))

EXEDIR:=../$(EXEDIR)
MAKEOPT = -f ../Makefile
all: yes_opt

no_opt:
	$(MAKE) $(MAKEOPT) -C obj targ
yes_opt:
	$(MAKE) $(MAKEOPT) -C obj_opt EXESUFFIX=.opt OPT="$(OPT_FFLAGS)" targ

targ: $(TARG)

# dependencies
addrot.o:     param.f c_sys.f          c_ran.f
andersen.o:   param.f c_sys.f          c_ran.f
angmom.o:     param.f
atomdiatom.o: param.f                  c_ran.f
bsstep.o:     param.f
carttojac.o:  param.f
checkfrus.o:  param.f c_sys.f
checkhop.o:   param.f c_sys.f          c_ran.f
cn.o:         param.f
decocheck.o:  param.f c_sys.f          c_ran.f
detwell.o:    param.f
derivs.o:     param.f c_sys.f
diamin.o:     param.f
diapot.o:     param.f c_sys.f
diatom.o:     param.f                  c_ran.f
dint.o:       param.f c_sys.f c_traj.f c_ran.f
driver.o:     param.f c_sys.f c_traj.f c_ran.f c_output.f
ejscale.o:    param.f
elecdeco.o:   param.f c_sys.f
ems.o:        param.f c_sys.f c_traj.f c_ran.f c_output.f
ewkb.o:       param.f
finalstate.o: param.f c_sys.f c_traj.f c_ran.f c_output.f
getdvec2.o:   param.f
getgrad.o:    param.f c_sys.f
getpem.o:     param.f c_sys.f
getrho.o:     param.f
getrel.o:     param.f c_sys.f
getrelbz2.o:  param.f c_sys.f
gettemp.o:    param.f
header.o:
hessems.o:    param.f c_sys.f c_traj.f c_ran.f c_output.f
honey.o:      param.f                          c_output.f
hop.o:        param.f c_sys.f          c_ran.f
initelec.o:   param.f c_sys.f c_traj.f
initmol.o:    param.f c_sys.f c_traj.f
lindemann.o:  param.f
mmid.o:       param.f
noang.o:      param.f
normod.o:     param.f
nmturn.o:     param.f
nmpot.o:      param.f
orient.o:     param.f c_sys.f c_traj.f c_ran.f
period.o:     param.f
pjsplit.o:    param.f
popnorm.o:    param.f c_sys.f          c_ran.f
preatomdiatom.o: param.f
premol.o:     param.f c_sys.f
pzextr.o:
radialdist.o: param.f                          c_output.f
rancol.o:     param.f c_sys.f          c_ran.f
ranclu.o:     param.f c_sys.f          c_ran.f
ranj.o:       param.f c_sys.f          c_ran.f
rantherm.o:   param.f c_sys.f          c_ran.f
ranrot.o:     param.f c_sys.f          c_ran.f
ransamp.o:    param.f c_sys.f          c_ran.f
ranstate.o:                            c_ran.f
ranwell.o:    param.f                  c_ran.f
rarray.o:     param.f
readin.o:     param.f c_sys.f c_traj.f c_ran.f
rk4.o:        param.f
rtox.o:       param.f
spin.o:       param.f                  c_ran.f
stodeco.o:    param.f c_sys.f          c_ran.f
takestep.o:   param.f c_sys.f c_traj.f
timing.o:
turn.o:       param.f
vibwells.o:   param.f c_sys.f c_traj.f c_ran.f
vwkb.o:       param.f
xptoy.o:      param.f

COMPILE = $(FC) $(FFLAGS) $(OPT)

$(TARG): $(OBJ)
	$(LD) $(LDFLAGS) $(OBJ) $(LIBS) -o $(TARG)

%.o : %.f
	$(COMPILE) $< -o $@

%.o : %.F
	$(COMPILE) $(FTPPFLAGS) $< -o $@

clean:
	-rm -f obj/* obj_opt/*

realclean:
	-rm -f obj/* obj_opt/* ../exe/*

