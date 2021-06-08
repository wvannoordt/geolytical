#ifndef GEO_SPHERE_H
#define GEO_SPHERE_H
#include "AnalyticalGeometry.h"
#include "DeformationTypes.h"
#include "Icosohedron.h"
namespace geolytical
{
    class Sphere : public Icosohedron
    {
        public:
            Sphere(int nxEquator_in, bbox bounds_in);
            ~Sphere(void);
        private:
            void MapPointsToSphere(void);
    };
}
#endif