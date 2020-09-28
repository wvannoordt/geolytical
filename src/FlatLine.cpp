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
        bounds = bounds_in;
        nx = nx_in+1;
        numPoints = nx+2;
        nFaces = numPoints;
        nSize = 3*nFaces;
        dealloc = true;
        points = (double*)malloc(numPoints*3*sizeof(double));
        CreatePoints();
    }
    
    void FlatLine::CreatePoints(void)
    {
        double dx = (bounds.xmax-bounds.xmin)/(nx-1);
        for (int i = 0; i < nx; i++)
        {
            *(points+3*i + 0) = bounds.xmin+i*dx;
            *(points+3*i + 1) = bounds.ymax;
            *(points+3*i + 2) = 0.0;
        }
        //LR
        *(points+3*nx +  0) = bounds.xmax;
        *(points+3*nx +  1) = bounds.ymin;
        *(points+3*nx +  2) = 0.0;
        //LL
        *(points+3*nx +  3) = bounds.xmin;
        *(points+3*nx +  4) = bounds.ymin;
        *(points+3*nx +  5) = 0.0;
    }
    
    FlatLine::~FlatLine(void)
    {
        if (dealloc)
        {
            dealloc = false;
            free(points);
        }
    }
    
    void FlatLine::OutputToVtk(std::string filename)
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
        for (int i = 0; i < nx-1; i++) myfile << "1" << std::endl;
        myfile << "1" << std::endl;
        myfile << "1" << std::endl;
        myfile << "1" << std::endl;
        myfile.close();
    }
}
