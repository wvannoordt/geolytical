#include <iostream>
#include "geolytical.h"
#include <cmath>
bbox bounds;

void defineBump(double* x, double* y)
{
    double x_prev = *x;
    if (x_prev*x_prev < 0.04)
    {
        *y = 0.05;
    }
}

int main(void)
{
    bounds.xmin = -1.0;
    bounds.xmax = 1.0;
    bounds.ymin = -0.1;
    bounds.ymax = 0.0;
    bounds.zmin = -1.0;
    bounds.zmax = 1.0;
    geolytical::FlatLine plate(100, bounds);
    Transformation2D bump(defineBump);
    plate.Deform(bump);
    plate.OutputToVtk("output.vtk");
    return 0;
}
