#ifndef FLATPLATE_H
#define FLATPLATE_H
#include "AnalyticalGeometry.h"
#include "DeformationTypes.h"
namespace geolytical
{
    class FlatPlate : public AnalyticalGeometry
    {
        public:
            FlatPlate(int nx_in, int ny_in, bbox bounds_in);
            ~FlatPlate(void);
            void OutputToVtk(std::string filename);
            void OutputToVtk(std::string filename, bool doScalarXYZ);
            void CreatePoints(void);
        private:
            int nx, nz;
            int nFaces, nSize;
            int* lookup;
    };
}
#endif