#include <iostream>
#include "geolytical.h"
#include <cmath>
bbox bounds;
double h, Lx, Ly, Lz;

void deform(double* x, double* y, double* z);
void deform2(double* x, double* y, double* z);
int main(void)
{
    h = 1.0;
    bounds.xmin = -1.0;
    bounds.xmax = 1.0;
    bounds.ymin = -0.1;
    bounds.ymax = 0.0;
    bounds.zmin = -1.0;
    bounds.zmax = 1.0;
    geolytical::FlatPlate plate(60, 60, bounds);
    geolytical::FlatPlate plate2(60, 60, bounds2);
    plate.Deform(deform);
    
    
    plate.OutputToVtk("output.vtk");
    return 0;
}

void deform(double* x, double* y, double* z)
{
    double xp = *x;
    double yp = *y;
    double zp = *z;
    double theta = 1.2*(bounds.xmax-xp)*(xp-bounds.xmin)*(bounds.zmax-zp)*(zp-bounds.zmin);
    double r = 0.5*xp*xp+0.5*zp*zp;
    double s = sin(theta);
    double c = cos(theta);
    *y = yp+theta;
    *x = c*xp+s*zp;
    *z = -s*xp+c*zp;
}

void deform2(double* x, double* y, double* z)
{
    double xp = *x;
    double yp = *y;
    double zp = *z;
    if (yp>0)
    {
        *x = xp + yp*sin(yp);
        *y = yp;
        *z = zp + yp*cos(yp);
    }
}
