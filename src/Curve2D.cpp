#include "Curve2D.h"

namespace geolytical
{
    Curve2D::Curve2D(double* data, int num)
    {
        dimension = 2;
        points = (double*)malloc(3*num*sizeof(double));
        numPoints = num;
        numFaces = numPoints;
        bounds.xmin =  1e30;
        bounds.ymin =  1e30;
        bounds.zmin =  1e30;
        bounds.xmax = -1e30;
        bounds.ymax = -1e30;
        bounds.zmax = -1e30;
        Allocate();
        for (int i = 0; i < num; i++)
        {
            AddPoint(data[2*i], data[2*i+1], 0.0);
            bounds.xmin = (bounds.xmin<(points[3*i+0]))?(bounds.xmin):(points[3*i+0]);
            bounds.ymin = (bounds.ymin<(points[3*i+1]))?(bounds.ymin):(points[3*i+1]);
            bounds.zmin = (bounds.zmin<(points[3*i+2]))?(bounds.zmin):(points[3*i+2]);
            bounds.xmax = (bounds.xmax>(points[3*i+0]))?(bounds.xmax):(points[3*i+0]);
            bounds.ymax = (bounds.ymax>(points[3*i+1]))?(bounds.ymax):(points[3*i+1]);
            bounds.zmax = (bounds.zmax>(points[3*i+2]))?(bounds.zmax):(points[3*i+2]);
        }
        
        for (int i = 0; i < numPoints-1; i++)
        {
            AddFace(i, i+1);
        }
        AddFace(numPoints-1, 0);
    }
    
    Curve2D::~Curve2D(void)
    {

    }
}