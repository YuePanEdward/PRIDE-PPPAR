# Makefile created by Yuanxin Pan (yxpan@whu.edu.cn) on Mar 9, 2018

all: get_ctrl

# ASSIGMENTS

FFLAGS = -g -O3
FTN = gfortran

LIB = ../lib/libcom.a  

get_ctrl : get_ctrl.f90 $(LIB)
	$(FTN) $(FFLAGS) get_ctrl.f90 $(LIB) -o get_ctrl
