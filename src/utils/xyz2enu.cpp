//
// xyz2enu.cpp
//
//    Copyright (C) 2018 by Wuhan University
//
//    This program is an open source software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License (version 3) as
//    published by the Free Software Foundation.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License (version 3) for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
// author: Yuanxin Pan
//
// Converting XYZ time-series in kin_file to ENU

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void help()
{
    printf("usage: xyz2enu kin_file [refx refy refz]\n"
           "output: mjd sod enu(m)\n"
           "note: avgx avgy avgz are used by default\n");
}

void avg(FILE *fp, const char *fmt, double avg[])
{
    char buf[256];
    double x, y ,z;
    double sumx=0, sumy=0, sumz=0;
    int count=0;
    while (fgets(buf,sizeof(buf),fp)!=NULL)
    {
        if (3 != sscanf(buf, fmt, &x, &y, &z))
            continue;
        ++count;
        sumx += x;
        sumy += y;
        sumz += z;
    }

    avg[0] = sumx/count;
    avg[1] = sumy/count;
    avg[2] = sumz/count;
}

class Ellipsoid
{
public:
    Ellipsoid() {}
    Ellipsoid(double _a, double _b): a(_a), b(_b) {
        e2 = ((a*a-b*b)/(a*a));
    }

    double a, b, e2;
};

static const Ellipsoid WGS84(6378137.0, 6356752.314245);

void xyz2blh(const double xyz[], double blh[])
{
    double x=xyz[0], y=xyz[1], z=xyz[2];
    double r = sqrt(x*x + y*y);
    blh[1] = atan2(y, x);

    double B=0.7, B0, H, N;
    double e2 = WGS84.e2;
    do {
        B0 = B;
        double sinB = sin(B);
        N = WGS84.a/sqrt(1 - e2*sinB*sinB);
        H = r/cos(B) - N;
        B = atan(z/r/(1 - e2*N/(N+H)));
    } while(fabs(B-B0) > 1.E-12);   // 1E-10 -> 0.64 mm 

    blh[0]=B;
    blh[2]=H;
}

void ecef2enu(const double vec[], const double blh_ref[], double enu[])
{
    double sinB = sin(blh_ref[0]);
    double cosB = cos(blh_ref[0]);
    double sinL = sin(blh_ref[1]);
    double cosL = cos(blh_ref[1]);

    double e = -sinL*vec[0] + cosL*vec[1];
    double n = -cosL*sinB*vec[0] - sinL*sinB*vec[1] + cosB*vec[2];
    double u =  cosL*cosB*vec[0] + sinL*cosB*vec[1] + sinB*vec[2];
    enu[0] = e; enu[1] = n; enu[2] = u;
}

// cvt "xyz (m)" to "enu (m)" 
int main(int argc, char *argv[])
{
    FILE *fp=NULL;
    const char *fmt = "%d %lf %lf %lf %lf";
    const char *fmt_avg = "%*d %*lf %lf %lf %lf";

    if (argc==2 || argc==5) {
        fp = fopen(argv[1], "r");
        if (fp==NULL) {
            fprintf(stderr, " error: no such file: %s \n", argv[1]);
            return 1;
        }
    } else {
        help();
        return 1;
    }

    double xyz_ref[3];
    if (argc == 2) {
        avg(fp, fmt_avg, xyz_ref);
        rewind(fp);
    }
    else {
        xyz_ref[0] = atof(argv[2]);
        xyz_ref[1] = atof(argv[3]);
        xyz_ref[2] = atof(argv[4]);
    }
    double blh_ref[3];
    xyz2blh(xyz_ref, blh_ref);

    char buf[256];
    int mjd;
    double sod;
    double xyz[3], dx[3], enu[3];
    while (fgets(buf, sizeof(buf), fp) != NULL)
    {
        if (5 != sscanf(buf, fmt, &mjd, &sod, xyz, xyz+1, xyz+2))
            continue;
        dx[0] = xyz[0] - xyz_ref[0];
        dx[1] = xyz[1] - xyz_ref[1];
        dx[2] = xyz[2] - xyz_ref[2];
        ecef2enu(dx, blh_ref, enu);
        fprintf(stdout, "%5d %8.2f %12.4f %12.4f %12.4f\n", mjd, sod, enu[0], enu[1], enu[2]);
    }

    fclose(fp);
    return 0;
}
