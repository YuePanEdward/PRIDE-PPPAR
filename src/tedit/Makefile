# Makefile created by Yuanxin Pan (yxpan@whu.edu.cn), on Mar 10, 2018

all: tedit

DIR_GUARD = @mkdir -p $(@D)
OBJ_PATH = .obj/

# ASSIGMENTS
FTN = gfortran 
FFLAGS = -g -O3 -fno-align-commons

LIBCOM = ../lib/libcom.a

SRC := $(wildcard *.f90)
OBJ := $(SRC:.f90=.o)
OBJ := $(addprefix $(OBJ_PATH), $(OBJ))

tedit : $(OBJ) $(LIBCOM)
	$(FTN) $(FFLAGS) $(OBJ) $(LIBCOM) -o tedit

$(OBJ_PATH)%.o : %.f90
	$(DIR_GUARD)
	$(FTN) $(FFLAGS) -c $< -o $@
	@echo -n "$(OBJ_PATH)$*.o : " > $(@:.o=.d)
	@sed -n "/include/s/include '\(.*\.h\)'/\1/p" $< | xargs >> $(@:.o=.d)

.PHONY : clean
clean :
	-rm tedit $(OBJ_PATH)*

-include $(OBJ:.o=.d)

