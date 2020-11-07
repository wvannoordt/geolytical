#include "ExtrudedCurve2D.h"

namespace geolytical
{
    ExtrudedCurve2D::ExtrudedCurve2D(int nz_in, double zmin, double zmax, double* data, int num)
    {
        dimension = 2;
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
        for (int iz = 0; iz < nz; iz++)
        {
            for (int i = 0; i < nr; i++)
            {
                points[iz*3*nr+3*i] = data[2*i];
                points[iz*3*nr+3*i+1] = data[2*i+1];
                points[iz*3*nr+3*i+2] = zmin+iz*dz;
                bounds.xmin = (bounds.xmin<(points[3*i+0]))?(bounds.xmin):(points[3*i+0]);
                bounds.ymin = (bounds.ymin<(points[3*i+1]))?(bounds.ymin):(points[3*i+1]);
                bounds.xmax = (bounds.xmax>(points[3*i+0]))?(bounds.xmax):(points[3*i+0]);
                bounds.ymax = (bounds.ymax>(points[3*i+1]))?(bounds.ymax):(points[3*i+1]);
            }
        }
        dealloc = true;
        numFaceEnd = nr-2;
        nFaces = 2*nr*nz_in + 2*numFaceEnd;
        nSize = nFaces*4;
    }
    
    void ExtrudedCurve2D::OutputToVtk(std::string filename)
    {
        size_t numFacesWritten = 0;
        std::ofstream myfile;
        myfile.open(filename.c_str());
        myfile << "# vtk DataFile Version 3.0" << std::endl;
        myfile << "vtk output" << std::endl;
        myfile << "ASCII" << std::endl;
        myfile << "DATASET POLYDATA" << std::endl;
        myfile << "POINTS " << numPoints << " float" << std::endl;
        for (int i = 0; i < numPoints; i++) myfile << points[3*i] << " " << points[3*i+1] << " " << points[3*i+2] << std::endl;
        myfile << "POLYGONS " << nFaces << " " << nSize << std::endl;
        for (int iz = 0; iz < nz-1; iz++)
        {
            for (int ir = 0; ir < nr-1; ir++)
            {
                myfile << "3 " << (nr*(iz) + (ir+1)) << " " << (nr*(iz) + (ir)) << " " << (nr*(iz+1) + (ir)) << std::endl;numFacesWritten++;
                myfile << "3 " << (nr*(iz) + (ir+1)) << " " << (nr*(iz+1) + (ir)) << " " << (nr*(iz+1) + (ir+1)) << std::endl;numFacesWritten++;
            }
            myfile << "3 " << (nr*(iz) + (0)) << " " << (nr*(iz) + (nr-1)) << " " << (nr*(iz+1) + (nr-1)) << std::endl;numFacesWritten++;
            myfile << "3 " << (nr*(iz) + (0)) << " " << (nr*(iz+1) + (nr-1)) << " " << (nr*(iz+1) + (0)) << std::endl;numFacesWritten++;
        }
        
        myfile << "3 " << (nr*(0) + (0)) << " " << (nr*(0) + (1)) << " " << (nr*(0) + (nr-1)) << std::endl;numFacesWritten++;
        myfile << "3 " << (nr*(nz-1) + (1)) << " " << (nr*(nz-1) + (0)) << " " << (nr*(nz-1) + (nr-1)) << std::endl;numFacesWritten++;
        
        int low = 1;
        int high = nr-1;
        while(low+1 < high-1)
        {
            myfile << "3 " << (nr*(0) + (low)) << " " << (nr*(0) + (low+1)) << " " << (nr*(0) + (high)) << std::endl;numFacesWritten++;
            myfile << "3 " << (nr*(0) + (high)) << " " << (nr*(0) + (low+1)) << " " << (nr*(0) + (high-1)) << std::endl;numFacesWritten++;
            myfile << "3 " << (nr*(nz-1) + (low+1)) << " " << (nr*(nz-1) + (low)) << " " << (nr*(nz-1) + (high)) << std::endl;numFacesWritten++;
            myfile << "3 " << (nr*(nz-1) + (high)) << " " << (nr*(nz-1) + (high-1)) << " " << (nr*(nz-1) + (low+1)) << std::endl;numFacesWritten++;
            high --;
            low ++;
        }
        if ((nr%2)==0)
        {
            myfile << "3 " << (nr*(0) + (low)) << " " << (nr*(0) + (nr/2)) << " " << (nr*(0) + (high)) << std::endl;numFacesWritten++;
            myfile << "3 " << (nr*(nz-1) + (nr/2)) << " " << (nr*(nz-1) + (low)) << " " << (nr*(nz-1) + (high)) << std::endl;numFacesWritten++;
        }
        myfile << "CELL_DATA " << nFaces << std::endl;
        myfile << "SCALARS Components int" << std::endl;
        myfile << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < nFaces; i++) myfile << "1" << std::endl;
        myfile.close();
    }
    
    ExtrudedCurve2D::~ExtrudedCurve2D(void)
    {
        if (dealloc)
        {
            dealloc = false;
            free(points);
        }
    }
}
