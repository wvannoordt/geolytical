#include "Sphere.h"

namespace geolytical
{
    Sphere::Sphere(int nxEquator_in, bbox bounds_in)
    {
        nFaceOnEdge = (int)(((double)(nxEquator_in))/5.0);
        numTrisOnFace = this->GetNumTrisOnFace(nFaceOnEdge);
        dimension = 3;
        CountFaces();
        CountPoints();
        Allocate();
        CreatePoints();
        CreateFaces();
        bounds.xmin = -1;
        bounds.xmax =  1;
        bounds.ymin = -1;
        bounds.ymax =  1;
        bounds.zmin = -1;
        bounds.zmax =  1;
        MapPointsToSphere();
        this->RemapBoundingBox(bounds_in);
    }
    
    void Sphere::MapPointsToSphere(void)
    {
        for (size_t i = 0; i < numPoints; i++)
        {
            double x = points[3*i+0];
            double y = points[3*i+1];
            double z = points[3*i+2];
            double mg = sqrt(x*x+y*y+z*z);
            points[3*i+0] /= mg;
            points[3*i+1] /= mg;
            points[3*i+2] /= mg;
        }
    }
    
    Sphere::~Sphere(void)
    {
        
    }
}