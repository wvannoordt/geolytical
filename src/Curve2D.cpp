#include "Curve2D.h"

namespace geolytical
{
    Curve2D::Curve2D(double* data, int num)
    {
        dimension = 2;
        points = (double*)malloc(3*num*sizeof(double));
        numPoints = num;
        bounds.xmin =  1e30;
        bounds.ymin =  1e30;
        bounds.zmin =  1e30;
        bounds.xmax = -1e30;
        bounds.ymax = -1e30;
        bounds.zmax = -1e30;
        for (int i = 0; i < num; i++)
        {
            points[3*i] = data[2*i];
            points[3*i+1] = data[2*i+1];
            points[3*i+2] = 0.0;
            bounds.xmin = (bounds.xmin<(points[3*i+0]))?(bounds.xmin):(points[3*i+0]);
            bounds.ymin = (bounds.ymin<(points[3*i+1]))?(bounds.ymin):(points[3*i+1]);
            bounds.zmin = (bounds.zmin<(points[3*i+2]))?(bounds.zmin):(points[3*i+2]);
            bounds.xmax = (bounds.xmax>(points[3*i+0]))?(bounds.xmax):(points[3*i+0]);
            bounds.ymax = (bounds.ymax>(points[3*i+1]))?(bounds.ymax):(points[3*i+1]);
            bounds.zmax = (bounds.zmax>(points[3*i+2]))?(bounds.zmax):(points[3*i+2]);
        }
        dealloc = true;
        nFaces = numPoints;
        nSize = nFaces*3;
    }
    
    void Curve2D::OutputToVtk(std::string filename)
    {
        std::ofstream myfile;
        myfile.open(filename.c_str());
        myfile << "# vtk DataFile Version 3.0" << std::endl;
        myfile << "vtk output" << std::endl;
        myfile << "ASCII" << std::endl;
        myfile << "DATASET POLYDATA" << std::endl;
        myfile << "POINTS " << numPoints << " float" << std::endl;
        for (int i = 0; i < numPoints; i++) myfile << points[3*i] << " " << points[3*i+1] << " " << points[3*i+2] << std::endl;
        myfile << "LINES " << nFaces << " " << nSize << std::endl;
        for (int i = 0; i < numPoints-1; i++)
        {
            myfile << "2 " << i << " " << i+1 << std::endl;
        }
        myfile << "2 " << numPoints-1 << " " << 0 << std::endl;
        myfile << "CELL_DATA " << nFaces << std::endl;
        myfile << "SCALARS Components int" << std::endl;
        myfile << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < numPoints-1; i++) myfile << "1" << std::endl;
        myfile << "1" << std::endl;
        myfile.close();
    }
    
    Curve2D::~Curve2D(void)
    {
        if (dealloc)
        {
            dealloc = false;
            free(points);
        }
    }
}