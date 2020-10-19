#ifndef ANAGEOMCYLINDER_H
#define ANAGEOMCYLINDER_H
#include "AnalyticalGeometry.h"
#include "DeformationTypes.h"
namespace geolytical
{
    class Cylinder : public AnalyticalGeometry
    {
        public:
            Cylinder(int nx_in_cyl, int nr_in_rad, int layersonFace, bbox bounds_in);
            ~Cylinder(void);
            void OutputToVtk(std::string filename);
            void OutputToVtk(std::string filename, bool doScalarXYZ);
            void CreatePoints(void);
            void CountPoints(void);
            void CreateFacePoints(void);
            void CreateCylPoints(void);
            void WritePoint(double x, double y, double z);
            void CountFaces(void);
            void WriteFace(bool reverseOrder, std::ofstream & myfile, int p1, int p2, int p3);
            void WriteFaceTriangles(std::ofstream & myfile);
            void WriteCylTriangles(std::ofstream & myfile);
            int GetLayerStart(int ilayer);
            int GetCylinderRingStart(int ilayer);
        private:
            int nx, nr, nz;
            int nFaces, nSize;
            int* lookup;
            int halvingFactor, layersonFace, pointsOnFace;
            int frequencyOfHalving;
            int* numberOfRadialPointsAtLevel;
            bool* isIrregularLayer;
            size_t pointIdx;
            int dummy;
    };
}
#endif