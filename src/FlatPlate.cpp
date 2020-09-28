#include "GeoTypes.h"
#include "FlatPlate.h"
#include "DeformationTypes.h"
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
        nFaces = 2*((nx-1)*(nz-1) + nx + nz + 1);
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
    void FlatPlate::OutputToVtk(std::string filename)
    {
        OutputToVtk(filename, false);
    }
    
    void FlatPlate::OutputToVtk(std::string filename, bool doScalarXYZ)
    {
        if (doScalarXYZ)
        {
            lookup = (int*)malloc(nFaces*sizeof(int));
        }
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
                myfile << "3 " << idx11 << " " << idx00 << " " << idx01 << std::endl;
                myfile << "3 " << idx10 << " " << idx00 << " " << idx11 << std::endl;
            }
        }
        int idx010 = 0;
        int idx011 = (nx-1)*nz + 0;
        int idx110 = 0*nz + (nz-1);
        int idx111 = (nx-1)*nz + (nz-1);
        int idx000 = numPoints - 4;
        int idx001 = numPoints - 1;
        int idx100 = numPoints - 3;
        int idx101 = numPoints - 2;
        
        for (int i = 0; i < nx-1; i++)
        {
            myfile << "3 " << (i)*nz + (0) << " " << (i+1)*nz + (0) << " " << idx000 <<  std::endl;
            myfile << "3 " << idx100 << " " << (i+1)*nz + (nz-1) << " " << (i)*nz + (nz-1) << std::endl;
        }
        myfile << "3 " << idx001 << " " << idx000 << " " << idx011 << std::endl;
        myfile << "3 " << idx111 << " " << idx100 << " " << idx101 << std::endl;
        
        for (int j = 0; j < nz-1; j++)
        {            
            myfile << "3 " << (0)*nz + (j) << " " << idx000 << " " << (0)*nz + (j+1) << std::endl;
            //myfile << "3 " << idx000 << " " << (0)*nz + (j) << " " << (0)*nz + (j+1) << std::endl;
            myfile << "3 " << (nx-1)*nz + (j+1) << " " << idx001 << " " <<  (nx-1)*nz + (j) << std::endl;
            //myfile << "3 " << idx001 << " " << (nx-1)*nz + (j+1) << " " << (nx-1)*nz + (j) << std::endl;
        }
        myfile << "3 " << idx001 << " " << idx111 << " " << idx101 << std::endl;
        myfile << "3 " << idx110 << " " << idx000 << " " << idx100 << std::endl;
        
        myfile << "3 " << idx000 << " " << idx101 << " " << idx100 << std::endl;
        myfile << "3 " << idx000 << " " << idx001 << " " << idx101 << std::endl;
        
        myfile << "CELL_DATA " << nFaces << std::endl;
        myfile << "SCALARS Components int" << std::endl;
        myfile << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < nFaces; i++) myfile << "1" << std::endl;
        
        if (doScalarXYZ)
        {
            std::cout << "not implemented" << std::endl;
            abort();
            //myfile << "CELL_DATA " << nFaces << std::endl;
            //myfile << "SCALARS x float" << std::endl;
            //myfile << "LOOKUP_TABLE default" << std::endl;
            for (int i = 0; i < nFaces; i++)
            {
                
            }
        }
        if (doScalarXYZ)
        {
            free(lookup);
        }
        myfile.close();
    }
}
