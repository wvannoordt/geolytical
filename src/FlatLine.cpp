#include "GeoTypes.h"
#include "FlatLine.h"
#include "DeformationTypes.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
namespace geolytical
{
    FlatLine::FlatLine(int nx_in, bbox bounds_in)
    {
        dimension = 2;
        bounds = bounds_in;
        nx = nx_in+1;
        numPoints = nx+2;
        numFaces = numPoints;
        Allocate();
        CreatePoints();
    }
    
    void FlatLine::CreatePoints(void)
    {
        double dx = (bounds.xmax-bounds.xmin)/(nx-1);
        for (int i = 0; i < nx; i++)
        {
            AddPoint(bounds.xmin+i*dx, bounds.ymax, 0.0);
        }
        AddPoint(bounds.xmax, bounds.ymin, 0.0);
        AddPoint(bounds.xmin, bounds.ymin, 0.0);
        for (int i = 0; i < numPoints-1; i++)
        {
            AddFace(i, i+1);
        }
        AddFace(numPoints-1, 0);
    }
    
    FlatLine::~FlatLine(void)
    {

    }
}
