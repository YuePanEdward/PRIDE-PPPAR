#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: Install PRIDE-PPPAR                                             ##
##                                                                           ##
##  AUTHOR : Yuanxin Pan    yxpan@whu.edu.cn                                 ##
##                                                                           ##
##  VERSION: ver 1.4                                                         ##
##                                                                           ##
##  DATE   : Dec-17, 2019                                                    ##
##                                                                           ##
##              @ GNSS RESEARCH CENTER, WUHAN UNIVERSITY, 2018               ##
##                                                                           ##
##    Copyright (C) 2018 by Wuhan University                                 ##
##                                                                           ##
##    This program is free software: you can redistribute it and/or modify   ##
##    it under the terms of the GNU General Public License (version 3) as    ##
##    published by the Free Software Foundation.                             ##
##                                                                           ##
##    This program is distributed in the hope that it will be useful,        ##
##    but WITHOUT ANY WARRANTY; without even the implied warranty of         ##
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          ##
##    GNU General Public License (version 3) for more details.               ##
##                                                                           ##
##    You should have received a copy of the GNU General Public License      ##
##    along with this program.  If not, see <https://www.gnu.org/licenses/>. ##
##                                                                           ##
###############################################################################

RED='\033[0;31m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check compiler
gfortran --version > /dev/null 2>&1
if [ $? -ne 0 ]; then
    printf "${RED}error:${NC} no compiler: gfortran\n"
    printf "${RED}error:${NC} PRIDE-PPPAR installation failed\n"; exit
fi
make --version > /dev/null 2>&1
if [ $? -ne 0 ]; then
    printf "${RED}error:${NC} GNU make not found\n"
    printf "${RED}error:${NC} PRIDE-PPPAR installation failed\n"; exit
fi

# Compilation & Installation
install_dir=${HOME}/.PRIDE_PPPAR_BIN
rm -rf "$install_dir"
cd src && make && make install \
    && cd .. \
    && mkdir -p $install_dir \
    && cp -f ./bin/* $install_dir \
    && cp -f ./src/lib/libpride_pppar.so $install_dir \
    && chmod 755 ./scripts/*.sh \
    && cp -f ./scripts/pride_pppar.sh $install_dir/pride_pppar \
    && cp -f ./scripts/rtk2xyz.sh ./scripts/leap.sh $install_dir
if [ $? -eq 0 ]; then
    grep "^export PATH=$install_dir:\$PATH" ${HOME}/.bashrc > /dev/null 2>&1
    [ $? -ne 0 ] && echo "export PATH=$install_dir:\$PATH" >> ${HOME}/.bashrc
    grep "^export LD_LIBRARY_PATH=$install_dir:\$LD_LIBRARY_PATH" ${HOME}/.bashrc > /dev/null 2>&1
    [ $? -ne 0 ] && echo "export LD_LIBRARY_PATH=$install_dir:\$LD_LIBRARY_PATH" >> ${HOME}/.bashrc
    grep "^export LD_LIBRARY_PATH=$install_dir:\$LD_LIBRARY_PATH" ${HOME}/.bash_profile > /dev/null 2>&1
    [ $? -ne 0 ] && echo "export LD_LIBRARY_PATH=$install_dir:\$LD_LIBRARY_PATH" >> ${HOME}/.bash_profile
fi

# Output
ls ${install_dir}/lsq > /dev/null 2>&1
if [ $? -eq 0 ]; then
    echo -e "\033[1;31m" && cat ./doc/logo && echo -e "$NC"
    printf "${BLUE}::${NC} PRIDE-PPPAR (v1.4) installation successfully completed!\n"
    printf "${BLUE}::${NC} executable binaries are copy to $install_dir\n"
    printf "${BLUE}::${NC} $install_dir added to PATH\n"
else
    printf "${RED}errror:${NC} PRIDE-PPPAR installation failed!\n"
    printf "${BLUE}Note:${NC} try replacing src/lib/libpride_pppar.so with .so file in src/lib/shared/\n"
    printf "      according to your OS and gfortran version\n"
fi
