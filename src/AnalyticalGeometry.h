#ifndef ANA_GEO_H
#define ANA_GEO_H
#include "GeoTypes.h"
#include "DeformationTypes.h"
#include <string>
namespace geolytical
{
    class AnalyticalGeometry
    {
        public:
            AnalyticalGeometry(void){}
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
        protected:
            int numPoints;
            double* points;
            bbox bounds;
            bool dealloc;
    };
}
#endif
