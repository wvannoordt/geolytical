#ifndef FLATLINE_H
#define FLATLINE_H
#include "AnalyticalGeometry.h"
#include "DeformationTypes.h"
namespace geolytical
{
    class FlatLine : public AnalyticalGeometry
    {
        public:
            FlatLine(int nx_in, bbox bounds_in);
            ~FlatLine(void);
            void CreatePoints(void);
        private:
            int nx;
    };
}
#endif