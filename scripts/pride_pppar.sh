#!/bin/bash

###############################################################################
##                                                                           ##
##  PURPOSE: GPS PPP&PPP-AR data processing with PRIDE-PPPAR                 ##
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


######################################################################
##                        Message Colors                            ##
######################################################################
NC='\033[0m'
RED='\033[0;31m'
BLUE='\033[1;34m'
CYAN='\033[0;36m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'

MSGERR="${RED}error:$NC"
MSGWAR="${YELLOW}warning:$NC"
MSGINF="${BLUE}::$NC"
MSGSTA="${BLUE}===>$NC"

USECACHE=NO # YES/NO (uppercase!)

######################################################################
##                     Funciton definations                         ##
######################################################################
main()
{
    CheckCmdArgs "$@" || exit 1
    CheckExecutables || exit 1

    local ctrl_file="$1"
    local ymd_start=(${2:0:4} ${2:4:2} ${2:6:2})
    local ymd_end=(${3:0:4} ${3:4:2} ${3:6:2})
    local AR="${4:0:1}"

    ctrl_file=$(readlink -f "$ctrl_file") # convert to absolute path

    # Output processing infomation
    echo -e "$MSGINF Processing date range: ${ymd_start[*]}  <==>  ${ymd_end[*]}"
    echo -e "$MSGINF Control file: ${ctrl_file}"
    echo -e "$MSGINF AR switch: ${AR}"

    # Processing day-by-day
    readonly local mjd_start=$(ymd2mjd ${ymd_start[*]})
    readonly local mjd_end=$(ymd2mjd ${ymd_end[*]})
    local work_dir=$(pwd) ydoy
    for mjd in $(seq $mjd_start $mjd_end)
    do
        cd "${work_dir}"
        ydoy=($(mjd2ydoy ${mjd}))
        mkdir -p ./${ydoy[0]}/${ydoy[1]}
        cd ./${ydoy[0]}/${ydoy[1]}
        if [ $? -eq 0 ]; then
            ProcessSingleDay $mjd "$ctrl_file" ${AR} || echo -e "$MSGWAR ${ydoy[*]} processing failed"
        else
            echo -e "$MSGERR no such directory: ${work_dir}/$year/$doy"
            echo -e "$MSGWAR skip processing: $year $doy"
        fi
    done
}

CheckCmdArgs() { # purpose: chech whether command line arguments are right
                 # usage  : CheckCmdArgs "$@"
    if [ $# -ne 4 ]; then
        PRIDE_PPPAR_Help
        return 1
    else
        local ctrl_file="$1"
        local ymd_start="$2"
        local ymd_end="$3"
        local AR="$4"
        if [ ! -e "$ctrl_file" ]; then
            echo -e "$MSGERR no control file: $ctrl_file"
            return 1
        elif [ ${#ymd_start} -ne 8 ]; then
            echo -e "$MSGERR bad start_date format: $ymd_start"
            echo -e "$MSGINF start_date format: YYYYMMDD"
            return 1
        elif [ ${#ymd_end} -ne 8 ]; then
            echo -e "$MSGERR bad end_date format: $ymd_start"
            echo -e "$MSGINF end_date format: YYYYMMDD"
            return 1
        elif [ ${AR,,} != y -a ${AR,,} != n ]; then
            echo -e "$MSGERR unknown AR option: $AR"
            echo -e "$MSGINF AR option: Y(y)/N(n)"
            return 1
        fi
    fi
}

CheckExecutables() { # purpose: check whether all needed executables are callable
                     # usage  : CheckExecutables
    echo -e "$MSGSTA CheckExecutables..."
    #local shell=$(ps -p $$ | awk 'NR==2 {print $4}')
    #echo $shell
    #if [ "$shell" != "bash" ]; then
    #    echo -e "$MSGERR pride_pppar should be executed with bash"
    #    return 1
    #fi
    if [ which lsq > /dev/null 2>&1 ]; then
        echo -e "$MSGERR lsq not found" && return 1
    fi
    if [ which arsig > /dev/null 2>&1 ]; then
        echo -e "$MSGERR arsig not found" && return 1
    fi
    if [ which tedit > /dev/null 2>&1 ]; then
        echo -e "$MSGERR tedit not found" && return 1
    fi
    if [ which redig > /dev/null 2>&1 ]; then
        echo -e "$MSGERR redig not found" && return 1
    fi
    if [ which get_ctrl > /dev/null 2>&1 ]; then
        echo -e "$MSGERR get_ctrl not found" && return 1
    fi
    if [ which mergesp3 > /dev/null 2>&1 ]; then
        echo -e "$MSGERR mergesp3 not found" && return 1
    fi
    if [ rnx2rtkp --help > /dev/null 2>&1 ]; then
        echo -e "$MSGERR rnx2rtkp not found" && return 1
    fi
    if [ awk --help > /dev/null 2>&1 ]; then
        echo -e "$MSGERR awk not found" && return 1
    fi
    if [ which wget > /dev/null 2>&1 ]; then
        echo -e "$MSGERR wget not found" && return 1
    fi
    if [ which readlink > /dev/null 2>&1 ]; then
        echo -e "$MSGERR readlink not found" && return 1
    fi
    echo -e "$MSGSTA CheckExecutables done"
}

PRIDE_PPPAR_Help() { # purpose: print usage for PRIDE_PPPAR
                     # usage  : PRIDE_PPPAR_Help
    echo " -----------------------------------------------------------------------"
    echo "  Purpose  :    GPS PPP&PPP-AR data processing with PRIDE-PPPAR"
    echo "  Usage    :    pride_pppar ctrl_file start_data end_date AR(Y/N)"
    echo "                   -- ctrl_file : configuration file of PRIDE-PPPAR"
    echo "                   -- start_date: start date for processing, format: YYYYMMDD"
    echo "                   -- end_date  : end date for procesing, format: YYYYMMDD"
    echo "                   -- AR(Y/N)   : Swithc of Ambiguity Resolution"
    echo "  Example  :    pride_pppar config 20160101 20160101 Y"
    echo "  Copyright:    GNSS Research Center, Wuhan University, 2018"
    echo " -----------------------------------------------------------------------"
}

ProcessSingleDay() { # purpose: process data of single day
                     # usage  : ProcessSingleDay mjd ctrl_file AR(y/n)
    local mjd=$1
    local ctrl_file="$2"
    local AR=$3

    local ydoy=($(mjd2ydoy $mjd))
    local year=${ydoy[0]}
    local doy=${ydoy[1]}
    local ymd=($(ydoy2ymd ${ydoy[*]}))
    local mon=${ymd[1]}
    local day=${ymd[2]}
    unset ydoy ymd

    echo -e "$MSGSTA ProcessSingleDay $year $doy..."

    # Copy a local ctrl_file
    cp -f "$ctrl_file" . || return 1
    ctrl_file=$(basename "$ctrl_file")

    # Create ctrl_file for current day
    sed -i -e "/^Session time/s/-YYYY-/$year/g" \
           -e "/^Session time/s/-MM-/$mon/g" \
           -e "/^Session time/s/-DD-/$day/g" "$ctrl_file"
    sed -i "/Rinex directory/s/-YEAR-/$year/g; s/-DOY-/$doy/g" "$ctrl_file"

    # Prepare products & tables
    local table_dir=$(get_ctrl "$ctrl_file" "Table directory")
    local product_dir=$(get_ctrl "$ctrl_file" "Sp3 directory")
    CopyTables "$table_dir" $mjd || return 1
    PrepareProducts $mjd "$product_dir" ${ctrl_file} ${AR}
    if [ $? -ne 0 ]; then
        echo -e "$MSGERR PrepareProducts failed"
        return 1
    fi

    # Download rinexnav file
    local rinex_dir=$(get_ctrl "$ctrl_file" "Rinex directory")
    local rinexnav="brdc${doy}0.${year:2:2}n"
    local urlnav="ftp://igs.gnsswhu.cn/pub/gps/data/daily/${year}/$doy/${year:2:2}n/${rinexnav}.Z"
    if [ ! -f "${rinex_dir}/${rinexnav}" ]; then
        echo -e "$MSGSTA Downloading ${rinexnav}..."
        WgetDownload "$urlnav"
        if [ $? -eq 0 ]; then
            uncompress -f ${rinexnav}.Z && mv ${rinexnav} ${rinex_dir}
            [ -f ${rinex_dir}/${rinexnav} ] && echo -e "$MSGSTA Downloading ${rinexnav} done" || return 1
        else
            echo -e "$MSGERR download rinexnav failed: ${rinex_dir}/${rinexnav}"
            return 1
        fi
    fi

    sites=($(awk '/^ [0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z] [KS]/ {print $1}' "$ctrl_file"))
    [ ${#sites[@]} -eq 0 ] && echo -e "$MSGWAR ${year} ${doy}: no site to be processed" && return 1
    for site in ${sites[*]}
    do
        ProcessSingleSite "$site" $year $doy "${ctrl_file}" $AR
        if [ $? -ne 0 ]; then
            echo -e "$MSGWAR ProcessSingleDay: skip processing $year $doy $site"
            local tp types=(kin pos rck ztd htg amb res stt con)
            for tp in ${types[*]}
            do
                rm -f ${tp}_${year}${doy}
            done
        fi
    done
    rm -f ${ctrl_file}
    echo -e "$MSGSTA ProcessSingleDay $year $doy done"
}

ProcessSingleSite() { # purpose: process data of single site
                      # usage  : ProcessSingleSite site year doy ctrl_file AR(Y/N)
    local site=$1
    local year=$2
    local doy=$3
    local ctrl_file=$4
    local AR=${5:0:1}

    echo -e "$MSGSTA ProcessSingleSite $year $doy $site..."

    # Data format conversion
    # ...
    # RINEX-OBS check
    local rinex_dir=$(get_ctrl "$ctrl_file" "Rinex directory")
    local rinexobs="${rinex_dir}/${site}${doy}0.${year:2:2}o"
    local rinexnav="${rinex_dir}/brdc${doy}0.${year:2:2}n"
    if [ ! -f "$rinexobs" ]; then
        echo -e "$MSGWAR ${rinexobs} doesn't exist" && return 1
    fi

    # Prepare initial site's position
    local xyz=($(sed -n "/^ $site/ {s/$site//p; q}" "./sit.xyz"))
    if [ ${#xyz[@]} -ne 3 ]; then
        echo -e "$MSGSTA Prepare initial position ${site}..."
        local initial_pos=($(ComputeInitialPos "$rinexobs" "$rinexnav"))
        if [ ${#initial_pos[@]} -ne 3 ]; then
            echo -e "$MSGERR ProcessSingleSite: no position: $site"
            return 1
        else
            xyz=(${initial_pos[*]})
            printf " %s%16.4f%16.4f%16.4f\n" ${site} ${initial_pos[*]} >> sit.xyz
            echo -e "$MSGSTA Prepare initial position ${site} done"
        fi
    fi

    echo -e "$MSGSTA Data pre-processing..."
    # Data preprocess
    local interval=$(get_ctrl "$ctrl_file" "Interval")
    local positon_mode=$(grep "^ $site [KS]" "$ctrl_file" | awk '{print $2}') # Static/Kinematic
    local cutoff_elev=$( grep "^ $site [KS]" "$ctrl_file" | awk '{print $5}') # int, degree
    local rhd_file="rhd_${year}${doy}_${site}"
    local ymd=($(ydoy2ymd $year $doy))
    local cmd=""
    if [ "$positon_mode" == S -o "$positon_mode" == K ]; then
        cmd="tedit ${rinexobs} -int ${interval} -rnxn ${rinexnav} -xyz ${xyz[*]} \
            -len 86400 -short 1200 -lc_check only -rhd ${rhd_file} -pc_check 300 \
            -elev ${cutoff_elev} -time ${ymd[*]} 0 0 0"
    else
        echo -e "$MSGERR ProcessSingleSite: unknown position mode: $site $positon_mode"
        return 1
    fi
    cmd=$(tr -s " " <<< "$cmd")
    ExecuteWithoutOutput "$cmd" || return 1
    echo -e "$MSGSTA Data pre-processing done"

    # Create site ctrl_file
    local tmp_ctrl=$(basename `mktemp -u`)
    tmp_ctrl=${tmp_ctrl/tmp/config}
    sed '/^\*NAME/ q' $ctrl_file > ${tmp_ctrl}
    grep "^ ${site} [KS]" "$ctrl_file" >> ${tmp_ctrl}
    echo "-Station used" >> ${tmp_ctrl}

    # Data clean (iteration)
    echo -e "$MSGSTA Data cleaning..."
    local short=$(echo $interval | awk '{printf("%.0f\n", 600/$1)}')
    local jumps=(400 200 100 50 50) # (500 300 200 100 50)
    for jump in ${jumps[*]}
    do
        cmd="lsq ${tmp_ctrl}"
        ExecuteWithoutOutput "$cmd" || return 1
        cmd="redig res_${year}${doy} -jmp $jump -sht $short"
        ExecuteWithoutOutput "$cmd" || return 1
    done
    # cmd="redig res_${year}${doy}"
    # Execute "$cmd" || return 1
    echo -e "$MSGSTA Data cleaning done"

    # Final process
    echo -e "$MSGSTA Final processing..."
    cmd="lsq ${tmp_ctrl}"
    Execute "$cmd" || return 1
    if [ $AR == Y -o $AR == y ]; then
        cmd="arsig ${tmp_ctrl}"
        Execute "$cmd" || return 1
        cmd="lsq ${tmp_ctrl}"
        Execute "$cmd" || return 1
    fi
    echo -e "$MSGSTA Final processing done"

    # Rename result files
    local tp types=(kin pos rck ztd htg amb res stt con) fn
    for tp in ${types[*]}
    do
        fn=${tp}_${year}${doy}
        [ -f ${fn} ] && mv -f ${fn} ${fn}_${site}
    done

    rm -f ${tmp_ctrl}
    echo -e "$MSGSTA ProcessSingleSite $year $doy $site done"
}

ComputeInitialPos() { # purpose: compute intial postion with rnx2rtkp
                      # usage  : ComputeInitialPos rinexobs rinexnav
    local rinexobs="$1"
    local rinexnav="$2"
    local site=$(basename "$rinexobs"); site=${site:0:4}
    local tmp_file=$(mktemp -u)

    rnx2rtkp -p 0 -ti 3600 -e -o ${tmp_file} ${rinexobs} ${rinexnav} || return 1
    [ ! -e ${tmp_file} ] && return 1
    awk -v nam=$site '{
        if (substr($0,1,1) != "%") {
             printf("%16.4f%16.4f%16.4f\n",$3,$4,$5);
             exit;
        }
    }' ${tmp_file} && rm ${tmp_file} || return 1
}

CopyTables() { # purpose: copy PRIDE-PPPAR needed tables to working directory
               # usage  : CopyTables table_dir mjd
    echo -e "$MSGSTA CopyTables..."
    local table_dir="$1"
    local mjd=$2
    local tables=(file_name jpleph_de405 leap.sec oceanload orography_ell sit.xyz)
    for table in ${tables[*]}
    do
        if [ ! -f "$table_dir/$table" ]; then
             echo -e "$MSGERR CopyTables: no such file: $table_dir/$table"
             return 1
        fi
        cp -f "$table_dir/$table" .
    done

    # IGS ATX
    local abs_atx=null
    if [ $mjd -lt 55668 ]; then
        abs_atx="igs05_1627.atx"
    elif [ $mjd -lt 57782 ]; then
        abs_atx="igs08_1930.atx"
    else
        abs_atx="igs14_2082.atx"
    fi
    cp "$table_dir/$abs_atx" ./abs_igs.atx
    [ $? -ne 0 ] && echo -e "$MSGERR CopyTables: no such file: $table_dir/$abs_atx" && return 1

    # Check valid time of table
    local mjd=`tail -2 leap.sec | head -1 | awk '{print $1}'`
    local ydoy=(`mjd2ydoy $((mjd))`)
    echo -e "$MSGINF leap.sec     is valid until ${ydoy[*]} (updated by leap.sh)"
    echo -e "$MSGINF jpleph_de405 is valid until 2040 007 (updated by PRIDELab website)"

    echo -e "$MSGSTA CopyTables done"
}

PrepareProducts() { # purpose: prepare PRIDE-PPPAR needed products in working directory
                    # usage  : PrepareProducts mjd products_dir AR(y/n)
    local mjd_mid=$1
    local products_dir="$2"
    local ctrl_file="$3"
    local AR=${4:0:1}
    [ -d $products_dir ] || mkdir -p "$products_dir"

    echo -e "$MSGSTA PrepareProducts..."

    local ydoy=($(mjd2ydoy $mjd_mid))
    local ymd=($(ydoy2ymd ${ydoy[*]}))
    local wkdow=($(mjd2wkdow $mjd_mid))
    local year=${ydoy[0]}
    local doy=${ydoy[1]}

    local clk="whp${wkdow[0]}${wkdow[1]}.clk.Z"
    [ $year -gt 2018 ] && clk="WHU5IGSFIN_${year}${ydoy[1]}0000_01D_30S_CLK.CLK.Z"
    local clk_url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/${year}/clock/${clk}"
    CopyOrDownloadProduct "$products_dir/$clk" "$clk_url" || return 1
    uncompress -f ${clk}

    local fcb="WHU0IGSFIN_${year}${ydoy[1]}0000_01D_01D_ABS.BIA.Z"
    local fcb_url="ftp://igs.gnsswhu.cn/pub/whu/phasebias/${year}/bias/${fcb}"
    CopyOrDownloadProduct "$products_dir/$fcb" "$fcb_url"
    if [ $? -ne 0 ]; then
        echo -e "$MSGWAR PrepareProducts: $fcb download failed"
        [ $AR == y -o $AR == Y ] && echo -e "$MSGERR no phase bias product: $fcb" && return 1
    else
        uncompress -f ${fcb}
    fi

    local lastmon=(`LastYearMonth ${ymd[*]}`)
    local dcb1="P1C1${year:2:2}${ymd[1]}_RINEX.DCB.Z"
    local dcb2="P2C2${year:2:2}${ymd[1]}_RINEX.DCB.Z"  # not necessary
    local dcb1url="ftp://ftp.aiub.unibe.ch/CODE/${year}/${dcb1}"
    local dcb2url="ftp://ftp.aiub.unibe.ch/CODE/${year}/${dcb2}"
    CopyOrDownloadProduct "$products_dir/$dcb1" "$dcb1url"
    if [ $? -ne 0 ]; then
        dcb1="P1C1${lastmon[0]:2:2}${lastmon[1]}_RINEX.DCB.Z"
        dcb2="P2C2${lastmon[0]:2:2}${lastmon[1]}_RINEX.DCB.Z"
        dcb1url="ftp://ftp.aiub.unibe.ch/CODE/${year}/${dcb1}"
        dcb2url="ftp://ftp.aiub.unibe.ch/CODE/${year}/${dcb2}"
        CopyOrDownloadProduct "$products_dir/$dcb1" "$dcb1url" || return 1
    fi

    uncompress -f ${dcb1}
    CopyOrDownloadProduct "$products_dir/$dcb2" "$dcb2url"
    if [ $? -ne 0 ]; then
        echo -e "$MSGWAR PrepareProducts: $dcb2 download failed"
    else
        uncompress -f ${dcb2}
    fi

    sed -n '8 p' ${fcb%.Z} | grep "rapid" > /dev/null 2>&1
    local rapid=$?  # whether use rapid products
    [ $rapid -eq 0 ] && echo -e "$MSGINF NOTE: Rapid Products Used"

    local erp erp_url
    if [ $rapid -ne 0 ]; then
        erp="COD${wkdow[0]}${wkdow[1]}.ERP.Z"
        erp_url="ftp://ftp.aiub.unibe.ch/CODE/${ydoy[0]}/${erp}"
        CopyOrDownloadProduct "$products_dir/$erp" "$erp_url"
        if [ $? -ne 0 ]; then
            erp="COD${wkdow[0]}7.ERP.Z"
            erp_url="ftp://ftp.aiub.unibe.ch/CODE/${ydoy[0]}/${erp}"
            CopyOrDownloadProduct "$products_dir/$erp" "$erp_url" || return 1
        fi
        uncompress -f ${erp}
    fi

    local sp3s erps tmpy
    local sp3 sp3_url i=0
    for mjd in $((mjd_mid-1)) mjd_mid $((mjd_mid+1))
    do
        tmpy=($(mjd2ydoy $mjd))
        wkdow=($(mjd2wkdow $mjd))

        if [ $rapid -eq 0 ]; then
            erp="COD${wkdow[0]}${wkdow[1]}.ERP_M.Z"
            erp_url="ftp://ftp.aiub.unibe.ch/CODE/${tmpy[0]}_M/$erp"
            CopyOrDownloadProduct "$products_dir/$erp" "$erp_url" || return 1
            uncompress -f ${erp}
            erps[$((i))]=${erp%.Z}

            sp3="COD${wkdow[0]}${wkdow[1]}.EPH_M.Z"
            sp3_url="ftp://ftp.aiub.unibe.ch/CODE/${tmpy[0]}_M/$sp3"
            CopyOrDownloadProduct "$products_dir/$sp3" "$sp3_url" || return 1
            uncompress -f ${sp3}
            sp3s[$((i++))]=${sp3%.Z}
        else
            sp3="COD${wkdow[0]}${wkdow[1]}.EPH.Z"
            sp3_url="ftp://ftp.aiub.unibe.ch/CODE/${tmpy[0]}/${sp3}"
            CopyOrDownloadProduct "$products_dir/$sp3" "$sp3_url" || return 1
            uncompress -f ${sp3}
            sp3s[$((i++))]=${sp3%.Z}
        fi
    done

    grep '^ [0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z][0-9a-zA-Z] .*VM1' ${ctrl_file} 2>&1 > /dev/null
    if [ $? -eq 0 ]; then
        echo -e "$MSGSTA Downloading VMF1 GRID..."
        local vmf vmf_url hour
        # Current Day (for interpolation)
        for hour in `seq 0 6 18 | awk '{printf("%02d\n",$1)}'`
        do
            vmf="VMFG_${ymd[0]}${ymd[1]}${ymd[2]}.H${hour}"
            vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/2.5x2/VMF1/VMF1_OP/${ydoy[0]}/${vmf}"
            CopyOrDownloadProduct "$products_dir/$vmf" "$vmf_url" || return 1
        done

        # Previous Day (for interpolation)
        tmpy=($(mjd2ydoy $((mjd_mid-1))))
        tmpy=($(ydoy2ymd ${tmpy[*]}))
        vmf="VMFG_${tmpy[0]}${tmpy[1]}${tmpy[2]}.H18"
        vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/2.5x2/VMF1/VMF1_OP/${tmpy[0]}/${vmf}"
        CopyOrDownloadProduct "$products_dir/$vmf" "$vmf_url" || return 1

        # Next Day (for interpolation)
        tmpy=($(mjd2ydoy $((mjd_mid+1))))
        tmpy=($(ydoy2ymd ${tmpy[*]}))
        vmf="VMFG_${tmpy[0]}${tmpy[1]}${tmpy[2]}.H00"
        vmf_url="http://vmf.geo.tuwien.ac.at/trop_products/GRID/2.5x2/VMF1/VMF1_OP/${tmpy[0]}/${vmf}"
        CopyOrDownloadProduct "$products_dir/$vmf" "$vmf_url" || return 1

        cat VMFG_* > vmf_${ydoy[0]}${ydoy[1]} || return 1
        echo -e "$MSGSTA Downloading VMF1 GRID done"
    fi

    # rename products
    mv ${clk%.Z} sck_${ydoy[0]}${ydoy[1]} || return 1
    [ -e ${fcb%.Z} ] && mv ${fcb%.Z} fcb_${ydoy[0]}${ydoy[1]}
    mv ${dcb1%.Z} P1C1.dcb           || return 1
    [ -e ${dcb2%.Z} ] && mv ${dcb2%.Z} P2C2.dcb

    # Generate igserp
    if [ $rapid -eq 0 ]; then
        cat ${erps[0]} > igserp && tail -1 ${erps[1]} >> \
            igserp && tail -1 ${erps[2]} >> igserp || return 1
    else
        mv ${erp%.Z} igserp || return 1
    fi

    # Generate binary sp3
    local cmd="mergesp3 ${sp3s[*]} orb_temp"
    ExecuteWithoutOutput "${cmd}" || return 1
    cmd="sp3orb orb_temp -cfg ${ctrl_file}"
    ExecuteWithoutOutput "${cmd}" && rm -f orb_temp || return 1

    echo -e "$MSGSTA PrepareProducts done"
}

CopyOrDownloadProduct() { # purpose: copy or download a product
                          # usage  : CopyOrDownloadProduct file url
    local file="$1"
    local url="$2"
    if [ ${USECACHE} = "YES" -a -f "$file" ]; then
        cp -f "$file" .
    else
        WgetDownload "$url" || return 1
        cp -f $(basename "$url") "$file"
        return 0
    fi
}

WgetDownload() { # purpose: download a file with wget
                 # usage  : WgetDownload url
    local url="$1"
    local args="-nv -nc -c -t 3 --connect-timeout=10 --read-timeout=60"
    cmd="wget ${args} ${url}"
    $cmd
    [ -e $(basename "${url}") ] && return 0 || return 1
}

LastYearMonth() { # purpose: get last year-month
                  # usage  : LastYearMonth year month
    local year=$1
    local mon=$((10#$2))
    [ $((mon-1)) -lt 1  ] && mon=12 && year=$((year-1)) || mon=$((mon-1))
    printf "%4d %02d\n" $year $mon
}

UncompressFile() { # purpose: uncompress a file automatically
                   # usage  : UncompressFile file del(Y/N)
    local file="$1"
    local del=$2
    # file $file
}

Execute() {
    local cmd="$1"
    # echo $cmd
    time=`date +'%Y-%m-%d %H:%M:%S'`
    $cmd
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}($time)${NC} ${CYAN}$cmd${NC} executed ok" # | tee -a $log
        return 0
    else
        echo -e "${RED}($time)${NC} ${CYAN}$cmd${NC} executed failed"  #| tee -a $log
        return 1
    fi
}

ExecuteWithoutOutput() {
    local cmd="$1"
    # echo $cmd
    time=`date +'%Y-%m-%d %H:%M:%S'`
    $cmd > /dev/null 2>&1
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}($time)${NC} ${CYAN}$cmd${NC} executed ok" # | tee -a $log
        return 0
    else
        echo -e "${RED}($time)${NC} ${CYAN}$cmd${NC} executed failed"  #| tee -a $log
        echo -e "$MSGINF Here is the output:\n"
        $cmd
        return 1
    fi
}


######################################################################
##                      Time Convert Funcitons                      ##
##  Author: Shuyin Mao      shuyinm@whu.edu.cn                      ##
######################################################################
ymd2mjd()
{
    local year=$1
    local mon=$((10#$2))
    local day=$((10#$3))
    [ $year -lt 100 ] && year=$((year+2000))
    if [ $mon -le 2 ];then
        mon=$(($mon+12))
        year=$(($year-1))
    fi
    local mjd=`echo $year | awk '{print $1*365.25-$1*365.25%1-679006}'`
    mjd=`echo $mjd $year $mon $day | awk '{print $1+int(30.6001*($3+1))+2-int($2/100)+int($2/400)+$4}'`
    #local mjd=$(bc <<< "$year*365.25 - $year*365.25 % 1 - 679006")
    #mjd=$(bc <<< "($mjd + (30.6001*($mon+1))/1 + 2 - $year/100 + $year/400 + $day)/1")
    echo $mjd
}

mjd2ydoy()
{
    local mjd=$1
    local year=$((($mjd + 678940)/365))
    local mjd0=$(ymd2mjd $year 1 1)
    local doy=$(($mjd-$mjd0))
    while [ $doy -le 0 ];do
        year=$(($year-1))
        mjd0=$(ymd2mjd $year 1 1)
        doy=$(($mjd-$mjd0+1))
    done
    printf "%d %03d\n" $year $doy
}

ymd2wkdow()
{
    local year=$1
    local mon=$2
    local day=$3
    local mjd0=44243
    local mjd=$(ymd2mjd $year $mon $day)
    local difmjd=$(($mjd-$mjd0-1))
    local week=$(($difmjd/7))
    local dow=$(($difmjd%7))
    echo $week $dow
}

mjd2wkdow()
{
    local mjd=$1
    local mjd0=44243
    local difmjd=$(($mjd-$mjd0-1))
    local week=$(($difmjd/7))
    local dow=$(($difmjd%7))
    echo $week $dow
}

ydoy2ymd()
{
    local iyear=$1
    local idoy=$((10#$2))
    local days_in_month=(31 28 31 30 31 30 31 31 30 31 30 31)
    local iday=0
    [ $iyear -lt 100 ] && iyear=$((iyear+2000))
    local tmp1=$(($iyear%4))
    local tmp2=$(($iyear%100))
    local tmp3=$(($iyear%400))
    if [ $tmp1 -eq 0 -a $tmp2 -ne 0 ] || [ $tmp3 -eq 0 ]; then
       days_in_month[1]=29
    fi
    local id=$idoy
    local imon=0
    local days
    for days in ${days_in_month[*]}
    do
        id=$(($id-$days))
        imon=$(($imon+1))
        if [ $id -gt 0 ]; then
            continue
        fi
        iday=$(($id + $days))
        break
    done
    printf "%d %02d %02d\n" $iyear $imon $iday
}


######################################################################
##                               Entry                              ##
######################################################################
main "$@"

