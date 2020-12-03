#ifndef CURVE2D_H
#define CURVE2D_H
#include "AnalyticalGeometry.h"
#include "DeformationTypes.h"
namespace geolytical
{
    class Curve2D : public AnalyticalGeometry
    {
        public:
            Curve2D(double* data, int num);
            ~Curve2D(void);
        private:
    };
}
#endif