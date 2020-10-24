#ifndef ANA_GEO_H
#define ANA_GEO_H
#include "GeoTypes.h"
#include "DeformationTypes.h"
#include <string>
#include <iostream>
#include <fstream>
namespace geolytical
{
    class AnalyticalGeometry
    {
        public:
            AnalyticalGeometry(void){dimension = -999;}
            ~AnalyticalGeometry(void){}
            virtual void OutputToVtk(std::string filename){}
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
            bbox bounds;
            bool dealloc;
            int dimension;
    };
}
#endif
