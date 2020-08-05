#ifndef FLATPLATE_H
#define FLATPLATE_H
#include "AnalyticalGeometry.h"
namespace geolytical
{
    class FlatPlate : public AnalyticalGeometry
    {
        public:
            FlatPlate(int nx_in, int ny_in, bbox bounds_in);
            ~FlatPlate(void);
            void OutputToVtk(std::string filename);
            void Deform(void (*deformer)(double*,double*,double*));
            void CreatePoints(void);
        private:
            int nx, nz;
            int nFaces, nSize;
    };
}
#endif