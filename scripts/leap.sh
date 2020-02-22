#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: Generate table/leap.sec from sopac leap.sec file                ##
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

wget --help > /dev/null 2>&1
if [ $? -ne 0 ]; then
	echo "error: wget not found"
	exit 1
fi

wget -O temp -nv ftp://garner.ucsd.edu/pub/gamit/tables/leap.sec
inp="temp"

echo "+leap sec" >  leap.sec
awk 'BEGIN { lp=21 } {
    	if(index($0,"!")!=0) {
    		printf("%6d%4d\n", int($1-2400001.0), lp);
    		lp = lp + 1
        }
    }' $inp >> leap.sec
echo "-leap sec" >> leap.sec

rm -f temp
