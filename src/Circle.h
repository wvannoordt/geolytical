#ifndef GEO_CIRCL_H
#define GEO_CIRCL_H

#include "AnalyticalGeometry.h"

namespace geolytical
{
    class Circle : public AnalyticalGeometry
    {
        public:
            Circle(int nr_in, bbox bounds_in);
        private:
            void CreatePoints(void);
            int nr;
    };
}

#endif