# Makefile created by Yuanxin Pan (yxpan@whu.edu.cn), on Nov 25, 2017

LIB := libcom.a
all : $(LIB)

DIR_GUARD = @mkdir -p $(@D)
OBJ_PATH = .obj/

# Flags
FTN = gfortran
FFLAGS = -g -O3

SRC := $(wildcard *.f90)
OBJ := $(SRC:.f90=.o)
OBJ := $(addprefix $(OBJ_PATH), $(OBJ))
 
$(LIB) : $(OBJ)
	ar -rv $(LIB) $(OBJ)

$(OBJ_PATH)%.o : %.f90
	$(DIR_GUARD)
	$(FTN) $(FFLAGS) -c $< -o $@
	@echo -n "$(OBJ_PATH)$*.o : " > $(@:.o=.d)
	@sed -n "/include/s/include '\(.*\.h\)'/\1/p" $< | xargs >> $(@:.o=.d)

.PHONY : clean
clean :
	-rm $(LIB) $(OBJ_PATH)*

-include $(OBJ:.o=.d)

