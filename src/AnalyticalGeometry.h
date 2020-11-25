#ifndef ANA_GEO_H
#define ANA_GEO_H
#include "GeoTypes.h"
#include "DeformationTypes.h"
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#ifndef BOUNDS_CHECK
#define BOUNDS_CHECK 0
#endif
// this is horrible, never write code like this.
#if(BOUNDS_CHECK)
#define CHECKP ErrIfOBPoint()
#define CHECKF ErrIfOBFace()
#define FULLP ErrIfNotAllPoint()
#define FULLF ErrIfNotAllFace()
#else
#define CHECKP
#define CHECKF
#define FULLP
#define FULLF
#endif
#define GEO_MIN(a,b) (((a)>(b))?(b):(a))
#define GEO_MAX(a,b) (((a)<(b))?(b):(a))
namespace geolytical
{
    class AnalyticalGeometry
    {
        public:
            AnalyticalGeometry(void)
            {
                dealloc = false;
                dimension = -999;
                numPoints=-1;
                numFaces=-1;
                fidx=0;
                pidx=0;
                doComponentId=true;
            }
            ~AnalyticalGeometry(void)
            {
                if (dealloc)
                {
                    free(points);
                    free(faces);
                    free(compID);
                }
            }
            virtual void Allocate(void)
            {
                dealloc = true;
                if (dimension<0) {std::cout << "dimension has not been set!" << std::endl; abort();}
                if (numPoints<0) {std::cout << "numPoints has not been set!" << std::endl; abort();}
                if (numFaces<0) {std::cout << "numFaces has not been set!" << std::endl; abort();}
                points = (double*)malloc(3*numPoints*sizeof(double));
                if (dimension==2) {faces = (int*)malloc(2*numFaces*sizeof(double));}
                else if (dimension==3) {faces = (int*)malloc(3*numFaces*sizeof(double));}
                else {std::cout << "dimension is neither 2 nor 3!" << std::endl; abort();}
                compID = (int*)malloc(numFaces*sizeof(int));
                for (int i = 0; i < numFaces; i++) compID[i] = 1;
            }
            virtual void AddPoint(double x, double y)
            {
                CHECKP;
                points[pidx] = x;
                pidx++;
                CHECKP;
                points[pidx] = y;
                pidx++;
            }
            virtual void AddPoint(double x, double y, double z)
            {
                CHECKP;
                points[pidx] = x;
                pidx++;
                CHECKP;
                points[pidx] = y;
                pidx++;
                CHECKP;
                points[pidx] = z;
                pidx++;
            }
            virtual void AddFace(int x, int y)
            {
                CHECKF;
                faces[fidx] = x;
                fidx++;
                CHECKF;
                faces[fidx] = y;
                fidx++;
            }
            virtual void AddFace(int x, int y, int z)
            {
                CHECKF;
                faces[fidx] = x;
                fidx++;
                CHECKF;
                faces[fidx] = y;
                fidx++;
                CHECKF;
                faces[fidx] = z;
                fidx++;
            }
            virtual void AddFaceSwitch(int x, int y, int z)
            {
                CHECKF;
                faces[fidx] = y;
                fidx++;
                CHECKF;
                faces[fidx] = x;
                fidx++;
                CHECKF;
                faces[fidx] = z;
                fidx++;
            }
            virtual void AddFace(int x, int y, int z, bool switchValues)
            {
                if (switchValues)
                {
                    AddFaceSwitch(x, y, z);
                }
                else
                {
                    AddFace(x, y, z);
                }
            }
            virtual void ErrIfOBPoint(void)
            {
                if (pidx>=3*numPoints) {std::cout << "point index out of bounds!" << std::endl; abort();}
            }
            virtual void ErrIfOBFace(void)
            {
                if (fidx>=dimension*numFaces) {std::cout << "face index out of bounds!" << std::endl; abort();}
            }
            virtual void ErrIfNotAllPoint(void)
            {
                if (pidx != 3*numPoints) {std::cout << "output vtk before all points filled (" << pidx << "/" << 3*numPoints << ")" << std::endl; abort();}
            }
            virtual void ErrIfNotAllFace(void)
            {
                if (fidx != dimension*numFaces) {std::cout << "output vtk before all faces filled (" << fidx << "/" << dimension*numFaces << ")" << std::endl; abort();}
            }
            virtual void OutputToVtk(std::string filename)
            {
                FULLP;
                FULLF;
                numFaces = (fidx/(dimension));
                std::ofstream myfile;
                myfile.open(filename.c_str());
                myfile << "# vtk DataFile Version 3.0" << std::endl;
                myfile << "vtk output" << std::endl;
                myfile << "ASCII" << std::endl;
                myfile << "DATASET POLYDATA" << std::endl;
                myfile << "POINTS " << numPoints << " float" << std::endl;
                for (int i = 0; i < numPoints; i++)
                {
                    myfile << points[3*i] << " " << points[3*i+1] << " " << points[3*i+2] << std::endl;
                }
                myfile << "POLYGONS " << numFaces << " " << (1+dimension)*numFaces << std::endl;
                for (int i = 0; i < numFaces; i++)
                {
                    myfile << dimension;
                    for (int d = 0; d < dimension; d++)
                    {
                        myfile << " " << faces[i*dimension+d];
                    }
                    myfile << std::endl;
                }
                if (doComponentId)
                {
                    myfile << "CELL_DATA " << numFaces << std::endl;
                    myfile << "SCALARS Components int" << std::endl;
                    myfile << "LOOKUP_TABLE default" << std::endl;
                    for (int i = 0; i < numFaces; i++) myfile << compID[i] << std::endl;
                }
                myfile.close();
            }
            virtual void Deform(Transformation3D deformer)
            {
                for (int i = 0; i < numPoints; i++)
                {
                    deformer(points+3*i, points+3*i+1, points+3*i+2);
                }
            }
            virtual void Deform(Transformation2D deformer)
            {
                for (int i = 0; i < numPoints; i++)
                {
                    deformer(points+3*i, points+3*i+1);
                }
            }
            virtual void OutputPointsAsCSV(std::string filename)
            {
                FULLP;
                std::ofstream myfile;
                myfile.open(filename.c_str());
                for (int i = 0; i < numPoints; i++)
                {
                    myfile << points[3*i] << ", " << points[3*i+1] << ", " << points[3*i+2] << std::endl;
                }
                myfile.close();
            }
        protected:
            int numPoints;
            double* points;
            int numFaces;
            int* faces;
            int* compID;
            bbox bounds;
            bool dealloc;
            int dimension;
            bool doComponentId;
            size_t fidx, pidx;
    };
}
#endif
