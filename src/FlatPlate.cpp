#include "GeoTypes.h"
#include "FlatPlate.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
namespace geolytical
{
    FlatPlate::FlatPlate(int nx_in, int nz_in, bbox bounds_in)
    {
        bounds = bounds_in;
        nx = nx_in+1;
        nz = nz_in+1;
        numPoints = nx*nz + 4;
        nFaces = 2*nx*nz+10;
        nSize = 4*nFaces;
        dealloc = true;
        points = (double*)malloc(numPoints*3*sizeof(double));
        CreatePoints();
    }
    
    void FlatPlate::CreatePoints(void)
    {
        double dx = (bounds.xmax-bounds.xmin)/(nx-1);
        double dy = (bounds.zmax-bounds.zmin)/(nz-1);
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < nz; j++)
            {
                int idx = i*nz + j;
                *(points+3*idx + 0) = bounds.xmin+i*dx;
                *(points+3*idx + 2) = bounds.zmin+j*dy;
                *(points+3*idx + 1) = bounds.ymax;
            }
        }
        //000
        *(points+3*nx*nz +  0) = bounds.xmin;
        *(points+3*nx*nz +  2) = bounds.zmin;
        *(points+3*nx*nz +  1) = bounds.ymin;
        //010
        *(points+3*nx*nz +  3) = bounds.xmin;
        *(points+3*nx*nz +  5) = bounds.zmax;
        *(points+3*nx*nz +  4) = bounds.ymin;
        //011
        *(points+3*nx*nz +  6) = bounds.xmax;
        *(points+3*nx*nz +  8) = bounds.zmax;
        *(points+3*nx*nz +  7) = bounds.ymin;
        //001
        *(points+3*nx*nz +  9) = bounds.xmax;
        *(points+3*nx*nz + 11) = bounds.zmin;
        *(points+3*nx*nz + 10) = bounds.ymin;
    }
    
    FlatPlate::~FlatPlate(void)
    {
        if (dealloc)
        {
            dealloc = false;
            free(points);
        }
    }
    
    void FlatPlate::Deform(void (*deformer)(double*,double*,double*))
    {
        for (int i = 0; i < numPoints; i++)
        {
            deformer(points+3*i, points+3*i+1, points+3*i+2);
        }
    }
    
    void FlatPlate::OutputToVtk(std::string filename)
    {
        std::ofstream myfile;
        myfile.open(filename.c_str());
        myfile << "# vtk DataFile Version 3.0" << std::endl;
        myfile << "vtk output" << std::endl;
        myfile << "ASCII" << std::endl;
        myfile << "DATASET POLYDATA" << std::endl;
        myfile << "POINTS " << numPoints << " float" << std::endl;
        for (int i = 0; i < numPoints; i++) myfile << points[3*i] << " " << points[3*i+1] << " " << points[3*i+2] << std::endl;
        myfile << "POLYGONS " << nFaces << " " << nSize << std::endl;
        for (int i = 0; i < nx-1; i++)
        {
            for (int j = 0; j < nz-1; j++)
            {
                int idx00 = i*nz + j;
                int idx01 = i*nz + (j+1);
                int idx10 = (i+1)*nz + j;
                int idx11 = (i+1)*nz + (j+1);
                myfile << "3 " << idx00 << " " << idx11 << " " << idx01 << std::endl;
                myfile << "3 " << idx00 << " " << idx10 << " " << idx11 << std::endl;
            }
        }
        int idx100 = 0;
        int idx101 = (nx-1)*nz + 0;
        int idx110 = 0*nz + (nz-1);
        int idx111 = (nx-1)*nz + (nz-1);
        int idx000 = numPoints - 4;
        int idx001 = numPoints - 1;
        int idx010 = numPoints - 3;
        int idx011 = numPoints - 2;
        //-z
        myfile << "3 " << idx000 << " " << idx011 << " " << idx001 << std::endl;
        myfile << "3 " << idx000 << " " << idx010 << " " << idx011 << std::endl;
        //-y
        myfile << "3 " << idx000 << " " << idx001 << " " << idx101 << std::endl;
        myfile << "3 " << idx000 << " " << idx101 << " " << idx100 << std::endl;
        //y
        myfile << "3 " << idx011 << " " << idx110 << " " << idx111 << std::endl;
        myfile << "3 " << idx011 << " " << idx010 << " " << idx110 << std::endl;
        //-x
        myfile << "3 " << idx000 << " " << idx100 << " " << idx110 << std::endl;
        myfile << "3 " << idx000 << " " << idx110 << " " << idx010 << std::endl;
        //x
        myfile << "3 " << idx111 << " " << idx101 << " " << idx001 << std::endl;
        myfile << "3 " << idx111 << " " << idx001 << " " << idx011 << std::endl;
        
        myfile << "CELL_DATA " << nFaces << std::endl;
        myfile << "SCALARS Components int" << std::endl;
        myfile << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < nFaces; i++) myfile << "1" << std::endl;
        myfile.close();
    }
}