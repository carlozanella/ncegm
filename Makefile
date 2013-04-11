SRCS_f90d1 = \
demo.f90 \
newton.f90 \
ncegm.f90 \
ndifferential.f90 \
normal_cdf.f90 \
main.f90 \
kinds.f90 \
tauchen.f90 \
grids.f90 

OBJS_f90d1 = \
demo.o \
newton.o \
ncegm.o \
ndifferential.o \
normal_cdf.o \
main.o \
kinds.o \
tauchen.o \
grids.o 

SRC_DIR_f90d1 = 
OBJS_DIR = obj/
EXE_DIR = bin/

EXE = dcegm
FC = gfortran
IDIR = 
CFLAGS = -std=f2003 -ffree-line-length-none -W -fexpensive-optimizations -Ofast  -J$(OBJS_DIR) $(IDIR)
LFLAGS = -s 
LIBS = 

VPATH = $(SRC_DIR_f90d1):$(OBJS_DIR)
OBJS = $(addprefix $(OBJS_DIR), $(OBJS_f90d1))

all : $(EXE)

$(EXE) : $(OBJS_f90d1)
	@mkdir -p $(EXE_DIR)
	$(FC) -o $(EXE_DIR)$(EXE) $(OBJS) $(LFLAGS) $(LIBS)

$(OBJS_f90d1):
	@mkdir -p $(OBJS_DIR)
	$(FC) $(CFLAGS) -c $(SRC_DIR_f90d1)$(@:.o=.f90) -o $(OBJS_DIR)$@

clean :
	rm -f $(OBJS_DIR)*.*
	rm -f $(EXE_DIR)$(EXE)

# Dependencies of files
demo.o: \
    demo.f90 \
    grids.o \
    kinds.o \
    ncegm.o \
    ndifferential.o \
    tauchen.o
newton.o: \
    newton.f90 \
    kinds.o
ncegm.o: \
    ncegm.f90 \
    kinds.o \
    ndifferential.o
ndifferential.o: \
    ndifferential.f90 \
    kinds.o
normal_cdf.o: \
    normal_cdf.f90 \
    kinds.o
main.o: \
    main.f90 \
    demo.o
kinds.o: \
    kinds.f90
tauchen.o: \
    tauchen.f90 \
    kinds.o \
    normal_cdf.o
grids.o: \
    grids.f90 \
    kinds.o

