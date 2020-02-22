#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: Generate table/sit.xyz with rnx2rtkp                            ##
##                                                                           ##
##  AUTHOR : Yuanxin Pan    yxpan@whu.edu.cn                                 ##
##                                                                           ##
##  VERSION: ver 1.4                                                         ##
##                                                                           ##
##  DATE   : Jul-16, 2019                                                    ##
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


######################################################################
##                        Message Colors                            ##
######################################################################
NC='\033[0m'
RED='\033[0;31m'
MSGERR="${RED}error:$NC"

######################################################################
##                               Entry                              ##
######################################################################
# Check command line args
if [ $# -ne 1 ]; then
    echo "usage: rtk2xyz.sh rnxobs_directory"
    exit 1
elif [ ! -d "$1" ]; then
    echo -e "$MSGERR no such directory: $1"
    exit 1
fi
# Check executable
rnx2rtkp --help > /dev/null 2>&1
if [ $? -ne 0 ]; then
	echo -e "$MSGERR rnx2rtkp not found"
	exit 1
fi
wget --help > /dev/null 2>&1
if [ $? -ne 0 ]; then
	echo -e "$MSGERR wget not found"
	exit 1
fi

wk_dir=$(pwd)
cd $1

# Sites
sit=($(ls -1 *o | cut -c1-4))
if [ ${#sit[@]} -eq 0 ]; then
	echo "no rnxobs files"
    exit 0
fi

# Check if broadcast ephemeris exist
tmp=$(ls -1 *o | awk '{print $0;exit}')
doy=$(echo $tmp | cut -c5-7)
yr=$(echo $tmp | cut -c10-11)
if [ -e brdc${doy}0.${yr}n.Z ]; then
    gunzip brdc${doy}0.${yr}n.Z
fi
if [ ! -e brdc${doy}0.${yr}n ]; then
    wget -nv -nc -t 3 --connect-timeout=10 --read-timeout=60 ftp://cddis.gsfc.nasa.gov/pub/gps/data/daily/20$yr/$doy/${yr}n/brdc${doy}0.${yr}n.Z
    if [ "$?" -ne "0" ]; then
        echo -e "$MSGERR download rnxnav failed"
        exit 1
    fi
    gunzip brdc${doy}0.${yr}n.Z
fi

# Calculate site-by-site
rm -f sit.xyz sit.blh
for ss in ${sit[@]}; do
    echo $ss
    rnx2rtkp -o outxyz -p 0 -ti 3600 -e ${ss}${doy}0.${yr}o brdc${doy}0.${yr}n
    rnx2rtkp -o outblh -p 0 -ti 3600    ${ss}${doy}0.${yr}o brdc${doy}0.${yr}n
    # collect coordinates
    awk -v nam=$ss '{if(substr($0,1,1)!="%"){printf(" %16.6f%16.6f%10s\n",$4,$3,nam);exit}}' outblh >> sit.blh
    awk -v nam=$ss '{if(substr($0,1,1)!="%"){printf(" %s%16.4f%16.4f%16.4f\n",nam,$3,$4,$5);exit}}' outxyz >> sit.xyz
done
rm -f outxyz outblh
cd "$wk_dir"
