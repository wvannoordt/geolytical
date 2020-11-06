#ifndef EXTRUDEDCURVE2D_H
#define EXTRUDEDCURVE2D_H
#include "AnalyticalGeometry.h"
#include "DeformationTypes.h"
namespace geolytical
{
    class ExtrudedCurve2D : public AnalyticalGeometry
    {
        public:
            ExtrudedCurve2D(int nz_in, double zmin, double zmax, double* data, int num);
            ~ExtrudedCurve2D(void);
            void OutputToVtk(std::string filename);
        private:
            int nz, nr;
            int nFaces, nSize, numFaceEnd;
    };
}
#endif