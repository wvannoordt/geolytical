#include <iostream>
#include "geolytical.h"
#include <cmath>
bbox bounds;
double h, Lx, Ly, Lz;
#define XBL -0.62
#define XSTART 0
#define XEND 1
#define LDEV 0.9
#define LOUT 0.9
#define YPLATE 0.22
#define YMAX 1.16
#define ZMAX 0.128
#define NX 100
#define NZ 40

void deform(double* x, double* y, double* z);
int main(void)
{
    h = 1.0;
    bounds.xmin = XBL-LDEV;
    bounds.xmax = XEND+LOUT;
    bounds.ymin = -0.1;
    bounds.ymax = 0.0;
    bounds.zmin = 0.0;
    bounds.zmax = ZMAX;
    geolytical::FlatPlate plate(512, 256, bounds);
    plate.AddDoubleScalar("x", [](double x, double y, double z) {return x;});
    plate.AddDoubleScalar("y", [](double x, double y, double z) {return y;});
    plate.AddDoubleScalar("z", [](double x, double y, double z) {return z;});
    plate.Deform(deform);
    plate.OutputToVtk("output.vtk");
    return 0;
}

void deform(double* x, double* y, double* z)
{
    double xp = *x;
    double yp = *y;
    double zp = *z;
    if (xp < XSTART && yp > -0.008)
    {
        *y = YPLATE;
    }
    if (xp >= XSTART && yp > -0.008 && xp <= XEND)
    {
        *y = YPLATE*(1.0 - 10*xp*xp*xp + 15*xp*xp*xp*xp - 6*xp*xp*xp*xp*xp);
    }
}
