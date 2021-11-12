#include <iostream>
#include "geolytical.h"
#include <cmath>
using geolytical::v3d;

int main(void)
{
	bbox bounds;
    bounds.xmin = -1.0;
    bounds.xmax = 1.0;
    bounds.ymin = -0.1;
    bounds.ymax = 0.0;
    bounds.zmin = -1.0;
    bounds.zmax = 1.0;

    int nx = 100;
    geolytical::FlatPlate top(nx, nx, bounds);
    geolytical::FlatPlate bottom(nx, nx, bounds);
    
    v3d<> axis(0, 0, 1);
    v3d<> point(0.5*(bounds.xmin+bounds.xmax), bounds.ymax, 0.5*(bounds.zmin+bounds.zmax));
    const double pi = 3.14159265359;
    double theta = pi;
    
    top.Rotate(axis, point, theta);
    
    double channelHeight = 1.0;
    
    auto translate = [&](double* x, double* y, double* z) -> void {*y += channelHeight;};
    top.Transform(translate);
    
    bottom += top;
    
    bottom.OutputToVtk("output.vtk");
    return 0;
}
