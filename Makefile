#-----------------------------------------------------
# Makefile to compile the TGLFEP_driver system.
#-----------------------------------------------------

# Define compilers and flags.
# Compilers and flags

include ${GACODE_ROOT}/shared/install/make.inc.${GACODE_PLATFORM}
export EXTRA_LIBS = \
        ${GACODE_ROOT}/tglf/src/tglf_lib.a \
        ${GACODE_ROOT}/shared/harvest_client/libharvest.a \
        ${GACODE_ROOT}/shared/math/math_lib.a \
        ${GACODE_ROOT}/EPtran/EPtran_lib.a

ifeq ($(OPT),debug)
   FFLAGS=${FDEBUG}
else
   FFLAGS=${FOPT}
endif

EXEC=TGLFEP_driver

LLIB=TGLFEP_lib

OBJECTS = TGLFEP_interface.o\
          EPstd_tglf_map.o\
          EPGYRO_tglf_map.o\
          EPtran_tglf_map.o\
          TGLFEP_tglf_map.o\
          TGLFEP_ky.o\
          TGLFEP_TM.o\
          TGLFEP_ky_widthscan.o\
          TGLFEP_nEPscan.o\
          TGLFEP_ky_nEPscan.o\
          TGLFEP_bisection.o\
          TGLFEP_scalefactor.o\
          TGLFEP_mainsub.o 

.SUFFIXES : .o .f90 .f

all: $(LLIB).a $(EXEC)

$(EXEC): $(LLIB).a $(EXEC).o ${EXTRA_LIBS}
	$(FC) $(FFLAGS) -o $(EXEC) $(EXEC).o $(LLIB).a $(EXTRA_LIBS) ${LMATH}

$(LLIB).a: $(OBJECTS)
	$(ARCH) $(LLIB).a $(OBJECTS)

.f90.o :
	$(FC) $(FFLAGS) $(FMATH) -c $<

.f.o :
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.mod $(EXEC) $(LLIB).a
