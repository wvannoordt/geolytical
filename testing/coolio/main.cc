#include "geolytical.h"
void deformer(double* x, double* y, double* z)
{
    double xprev = *x;
    double yprev = *y;
    double zprev = *z;
    double r = exp(-0.04*xprev*xprev);
    double deltaz = r*sin(xprev);
    double deltay = r*cos(xprev);
    *y = yprev + deltay;
    *z = zprev + deltaz;
}
double divDelta(double x, double y, double z)
{
    double xprev = x;
    double div = 3.0 + exp(-0.04*xprev*xprev)*cos(xprev) + sin(xprev)*exp(-0.04*xprev*xprev)*(-0.08*xprev) - exp(-0.04*xprev*xprev)*sin(xprev) + cos(xprev)*exp(-0.04*xprev*xprev)*(-0.08*xprev);
    return div*div*cos(3*z)*sin(2*y);
}
int main(void)
{
    bbox bounds;
    bounds.xmin = -8.0;
    bounds.xmax = 8.0;
    bounds.ymin = 0.0;
    bounds.ymax = 2.0;
    bounds.zmin = 0.0;
    bounds.zmax = 2.0;
    geolytical::Cylinder cyl(120, 102, 46, bounds);
    cyl.Deform(deformer);
    cyl.AddDoubleScalar("coolVar", divDelta);
    cyl.OutputToVtk("output.vtk");
    return 0;
}
