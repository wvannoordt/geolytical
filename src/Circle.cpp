#include "Circle.h"
#include <cmath>
namespace geolytical
{
    Circle::Circle(int nr_in, bbox bounds_in)
    {
        dimension = 2;
        bounds = bounds_in;
        nr = nr_in;
        numPoints = nr;
        numFaces = numPoints;
        Allocate();
        CreatePoints();
    }
    
    void Circle::CreatePoints(void)
    {
        double dtheta = 2.0*pi/nr;
        double x0 = 0.5*(bounds.xmin + bounds.xmax);
        double y0 = 0.5*(bounds.ymin + bounds.ymax);
        for (int i = 0; i < nr; i++)
        {
            double x = cos(-i*dtheta);
            double y = sin(-i*dtheta);
            double rx = 0.5*(bounds.xmax - bounds.xmin);
            double ry = 0.5*(bounds.ymax - bounds.ymin);
            AddPoint(x0+rx*x, y0+ry*y, 0.0);
        }
        for (int i = 0; i < numPoints-1; i++)
        {
            AddFace(i, i+1);
        }
        AddFace(numPoints-1, 0);
    }
}