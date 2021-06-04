#ifndef ANA_GEO_H
#define ANA_GEO_H
#include "GeoTypes.h"
#include "DeformationTypes.h"
#include <string>
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <map>
#include "SurfaceVar.h"
#ifndef BOUNDS_CHECK
#define BOUNDS_CHECK 0
#endif
// this is horrible, never write code like this.
#if(BOUNDS_CHECK)
#define CHECKP ErrIfOBPoint()
#define CHECKF ErrIfOBFace()
#define FULLP ErrIfNotAllPoint()
#define FULLF ErrIfNotAllFace()
#else
#define CHECKP
#define CHECKF
#define FULLP
#define FULLF
#endif
#define GEO_MIN(a,b) (((a)>(b))?(b):(a))
#define GEO_MAX(a,b) (((a)<(b))?(b):(a))
namespace geolytical
{
    class AnalyticalGeometry
    {
        public:
            AnalyticalGeometry(void);
            ~AnalyticalGeometry(void);
            bool HasScalar(std::string name);
            virtual void Allocate(void);
            virtual void AddPoint(double x, double y);
            virtual void AddPoint(double x, double y, double z);
            virtual void AddFace(int x, int y);
            virtual void AddFace(int x, int y, int z);
            virtual void AddFaceSwitch(int x, int y, int z);
            virtual void AddFace(int x, int y, int z, bool switchValues);
            virtual void ErrIfOBPoint(void);
            virtual void ErrIfOBFace(void);
            virtual void ErrIfNotAllPoint(void);
            virtual void ErrIfNotAllFace(void);
            void RemoveScalar(std::string name);
            void GetCentroid(int faceIdx, double* xout, double* yout, double* zout);
            virtual SurfaceVar& AddDoubleScalar(std::string name, double value);
            virtual SurfaceVar& AddDoubleScalar(std::string name, double(*func)(int));
            virtual SurfaceVar& AddDoubleScalar(std::string name, double(*func)(double, double, double));
            virtual SurfaceVar& AddIntegerScalar(std::string name, int value);
            virtual SurfaceVar& AddIntegerScalar(std::string name, int (*func)(int));
            virtual SurfaceVar& AddIntegerScalar(std::string name, int (*func)(double, double, double));
            virtual void OutputToVtk(std::string filename);
            virtual void Deform(Transformation3D deformer);
            virtual void Deform(Transformation2D deformer);
            virtual void OutputPointsAsCSV(std::string filename);
            void ResetFaceCounter(void) {fidx = 0;}
        protected:
            int numPoints;
            double* points;
            int numFaces;
            int* faces;
            int* compID;
            bbox bounds;
            bool dealloc;
            int dimension;
            bool doComponentId;
            bool hasAnyScalars;
            size_t fidx, pidx;
            std::map<std::string, double*> doubleScalars;
            std::map<std::string, int*> integerScalars;
            std::map<std::string, SurfaceVar*> variableObjects;
    };
}
#endif
