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
            void CountPoints(void);
            void CountFaces(void);
            void CreatePoints(void);
            void CreateFaces(void);
            void GetLineArray(bool isX, bool upperFace, int level);
        private:
            int nx, nz;
            int* levelx;
            int* levelz;
            int* levelPtStart;
            int* numPointsInLayer;
            int* minPointIdArray;
            int numLevels;
            int* upperLineArray;
            int nLower;
            int* lowerLineArray;
            int nUpper;
            bool deleteLevels;
    };
}
#endif