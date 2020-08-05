#ifndef ANA_GEO_H
#define ANA_GEO_H
#include "GeoTypes.h"
#include <string>
namespace geolytical
{
    class AnalyticalGeometry
    {
        public:
            AnalyticalGeometry(void){}
            ~AnalyticalGeometry(void){}
            virtual void OutputToVtk(std::string filename){}
            virtual void Deform(void (*deformer)(double*,double*,double*)){}
        protected:
            int numPoints;
            double* points;
            bbox bounds;
            bool dealloc;
    };
}
#endif