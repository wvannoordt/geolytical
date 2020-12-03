#include "ExtrudedCurve2D.h"

namespace geolytical
{
    ExtrudedCurve2D::ExtrudedCurve2D(int nz_in, double zmin, double zmax, double* data, int num)
    {
        dimension = 3;
        nr = num;
        nz = nz_in+1;
        points = (double*)malloc(3*nr*nz*sizeof(double));
        numPoints = nz*num;
        bounds.xmin =  1e30;
        bounds.ymin =  1e30;
        bounds.zmin = zmin;
        bounds.xmax = -1e30;
        bounds.ymax = -1e30;
        bounds.zmax = zmax;
        double dz = (zmax-zmin)/nz_in;
        numFaceEnd = nr-2;
        numFaces = 2*nr*nz_in + 2*numFaceEnd;
        numPoints = nz * nr;
        Allocate();
        for (int iz = 0; iz < nz; iz++)
        {
            for (int i = 0; i < nr; i++)
            {
                double x = data[2*i];
                double y = data[2*i+1];
                double z = zmin+iz*dz;
                AddPoint(x, y, z);
                bounds.xmin = (bounds.xmin<(points[3*i+0]))?(bounds.xmin):(points[3*i+0]);
                bounds.ymin = (bounds.ymin<(points[3*i+1]))?(bounds.ymin):(points[3*i+1]);
                bounds.xmax = (bounds.xmax>(points[3*i+0]))?(bounds.xmax):(points[3*i+0]);
                bounds.ymax = (bounds.ymax>(points[3*i+1]))?(bounds.ymax):(points[3*i+1]);
            }
        }
        CreateFaces();
        
    }
    
    void ExtrudedCurve2D::CreateFaces(void)
    {
        for (int iz = 0; iz < nz-1; iz++)
        {
            for (int ir = 0; ir < nr-1; ir++)
            {
                AddFace((nr*(iz) + (ir+1)),(nr*(iz) + (ir)),(nr*(iz+1) + (ir)));
                AddFace((nr*(iz) + (ir+1)),(nr*(iz+1) + (ir)),(nr*(iz+1) + (ir+1)));
            }
            AddFace((nr*(iz) + (0)),(nr*(iz) + (nr-1)),(nr*(iz+1) + (nr-1)));
            AddFace((nr*(iz) + (0)),(nr*(iz+1) + (nr-1)),(nr*(iz+1) + (0)));
        }
        AddFace((nr*(0) + (0)), (nr*(0) + (1)), (nr*(0) + (nr-1)));
        AddFace((nr*(nz-1) + (1)), (nr*(nz-1) + (0)), (nr*(nz-1) + (nr-1)));
        
        int low = 1;
        int high = nr-1;
        while(low+1 < high-1)
        {
            AddFace((nr*(0) + (low))     , (nr*(0) + (low+1))    ,(nr*(0) + (high)));
            AddFace((nr*(0) + (high))    , (nr*(0) + (low+1))    ,(nr*(0) + (high-1)));
            AddFace((nr*(nz-1) + (low+1)), (nr*(nz-1) + (low))   ,(nr*(nz-1) + (high)));
            AddFace((nr*(nz-1) + (high)) , (nr*(nz-1) + (high-1)),(nr*(nz-1) + (low+1)));
            high --;
            low ++;
        }
        if ((nr%2)==0)
        {
            AddFace((nr*(0) + (low)),(nr*(0) + (nr/2)),(nr*(0) + (high)));
            AddFace((nr*(nz-1) + (nr/2)),(nr*(nz-1) + (low)),(nr*(nz-1) + (high)));
        }
    }
    
    ExtrudedCurve2D::~ExtrudedCurve2D(void)
    {

    }
}
