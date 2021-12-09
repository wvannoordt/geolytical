#include <iostream>
#include "geolytical.h"
#include <cmath>
#define NOSERADIUS (1.0)
#define AFTSCALE 5
#define XMAX 10
void noseSphere(double* x, double* y, double *z)
{
    double x_prev = *x;
    double y_prev = *y;
    double z_prev = *z;
    if (x_prev <= 1e-9)
    {
        double r2 = (y_prev*y_prev + z_prev*z_prev);
        if ((NOSERADIUS*NOSERADIUS)-r2 > 0)
        {
            *x = x_prev - sqrt((NOSERADIUS*NOSERADIUS)-r2);
        }
    }
}

void flare(double* x, double* y, double *z)
{
    double x_prev = *x;
    double y_prev = *y;
    double z_prev = *z;
    if (x_prev >= 1e-9)
    {
        *y = AFTSCALE*y_prev*x_prev/XMAX;
        *z = AFTSCALE*y_prev*x_prev/XMAX;
    }
}

int main(void)
{
	bbox bounds;
    bounds.xmin = -1.0;
    bounds.xmax = 1.0;
    bounds.ymin = -0.2;
    bounds.ymax = 0.2;
    bounds.zmin = -0.2;
    bounds.zmax = 0.2;
    geolytical::Cylinder cyl(120, 128, 46, bounds);

    // Transformation3D nose(noseSphere);
    // Transformation3D flareTransform(flare);
    // cyl.Deform(nose);
    //cyl.Deform(flareTransform);
    cyl.PermuteCoordinates(2, 1, 0);
    cyl.OutputToVtk("output.vtk");
    return 0;
}
