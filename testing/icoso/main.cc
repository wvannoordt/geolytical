#include <iostream>
#include "geolytical.h"
#include <cmath>

int main(void)
{
	bbox bounds;
    bounds.xmin = -0.47;
    bounds.xmax = 0.47;
    bounds.ymin = -0.47;
    bounds.ymax = 0.47;
    bounds.zmin = -0.47;
    bounds.zmax = 0.47;
    geolytical::Icosohedron ico(8, bounds);
    ico.OutputToVtk("sphere.vtk");
    return 0;
}
