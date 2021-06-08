#ifndef ICOSOHEDRON_H
#define ICOSOHEDRON_H
#include "AnalyticalGeometry.h"
#include "DeformationTypes.h"
namespace geolytical
{
    class Icosohedron : public AnalyticalGeometry
    {
        public:
            Icosohedron(int nxEquator_in, bbox bounds_in);
            Icosohedron(void){};
            ~Icosohedron(void);
            void CountPoints(void);
            void CountFaces(void);
            void CreatePoints(void);
            void CreateFaces(void);

        protected:
            int GetFacesOnEdge(int nxEquator_in);
            int GetNumTrisOnFace(int nFaceOnEdge_in);
            int nFaceOnEdge;
            int numTrisOnFace;
    };
}
#endif