#define PI 3.1415926535
#include <iostream>
#include "geolytical.h"
#include <cmath>
#include <string>
#include <cstdlib>
#include <vector>
int main(void)
{
    int nz = 162;
    int ndata_points = 258;
    double zmin = -0.05;
	double zmax = 0.25;
	double data [2*ndata_points];
#include "a_airfoil.hpp"
    double drmin = 5e-3;
    double* newPoints;
    int numNewPoints;
    int lasti=0;
    std::vector<int> idxes;
    for (int i = 0; i < ndata_points; i++)
    {
        double x1 = data[2*i];
        double x2 = data[2*lasti];
        double y1 = data[2*i+1];
        double y2 = data[2*lasti+1];
        double dr = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
        if (dr>drmin || i==0 || i==ndata_points/2)
        {
            lasti = i;
            idxes.push_back(i);
        }
    }
    
    newPoints = new double[2*idxes.size()];
    
    for (int i = 0; i < idxes.size(); i++)
    {
        newPoints[2*i] = data[2*idxes[i]];
        newPoints[2*i+1] = data[2*idxes[i]+1];
    }
    std::cout << "Reduced " << ndata_points-idxes.size() << " points" << std::endl;
    geolytical::ExtrudedCurve2D foil3D(nz, zmin, zmax, data, ndata_points);
    geolytical::Curve2D foil2D(data, ndata_points);
    foil3D.AddIntegerScalar("Components", [](double x, double y, double z){
    if (z < 0) return 1;
    if (z > 0.2) return 1;
    if (x < 0.1 && y > 0) return 2;
    if (x < 0.005) return 2;
    return 1;
    });
    foil2D.AddIntegerScalar("Components", [](double x, double y, double z){
    if (z < 0) return 1;
    if (z > 0.2) return 1;
    if (x < 0.1 && y > 0) return 2;
    if (x < 0.005) return 2;
    return 1;
    });
    foil2D.OutputToVtk("output2D.vtk");
    foil3D.OutputToVtk("output3D.vtk");
    
    delete[] newPoints;
    return 0;
}
