#define PI 3.1415926535
#include <iostream>
#include "geolytical.h"
#include <cmath>
#include <string>
#include <cstdlib>
int main(void)
{
    int nz = 162;
    double zmin = -0.05;
	double zmax = 0.25;
	double data [1024];
#include "data.hpp"
    geolytical::ExtrudedCurve2D foil3D(nz, zmin, zmax, data, 512);
	geolytical::Curve2D foil2D(data, 512);
	foil2D.OutputToVtk("output2D.vtk");
    foil3D.OutputToVtk("output3D.vtk");
	return 0;
}