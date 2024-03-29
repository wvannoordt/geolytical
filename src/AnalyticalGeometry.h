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
#include "v3d.h"
#include <type_traits>
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
    namespace PointCloudFormat
    {
        enum PointCloudFormat
        {
            csv,
            vtk
        };
    }
    
    class AnalyticalGeometry
    {
        public:
            AnalyticalGeometry(void);
            ~AnalyticalGeometry(void);
            bool HasScalar(std::string name);
            void Allocate(void);
            void AllocateFaces(void);
            void AllocatePoints(void);
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
            virtual SurfaceVar& GetScalar(std::string name);
            virtual SurfaceVar& GetVariable(std::string name);
            virtual void OutputToVtk(std::string filename);
            virtual void Deform(Transformation3D deformer);
            virtual void Deform(Transformation2D deformer);
            void OutputPointsAsCSV(std::string filename,  int numPointsToOutput, PointCloudFormat::PointCloudFormat format);
            void OutputPointsAsCSV(std::string filename,  PointCloudFormat::PointCloudFormat format);
            void OutputPointsAsCSV(std::string filename,  int numPointsToOutput);
            void OutputPointsAsCSV(std::string filename);
            void ResetFaceCounter(void) {fidx = 0;}
            void BufferFaces(void);
            void RemapBoundingBox(bbox newBox);
            virtual void PermuteCoordinates(int i1, int i2, int i3);
            double* GetPointBuffer(void) {return points;}
            int GetNumPoints(void) {return numPoints;}
            int GetNumFaces(void) {return numFaces;}
            v3d<double> GetCentroid(int faceIdx)
            {
                v3d<double> output(0.0,0.0,0.0);
                GetCentroid(faceIdx, &output[0], &output[1], &output[2]);
                return output;
            }
            template <class callable> void Transform(const callable& func)
            {
                for (int i = 0; i < numPoints; i++)
                {
                    func(points+3*i+0, points+3*i+1, points+3*i+2);
                }
            }
            
            template <class callable> void SetScalarValues(std::string variable, const callable& func)
            {
                
                SurfaceVar* var = NULL;
                bool isdouble = !std::is_integral<decltype(func(0,0,0))>::value;
                
                if (!this->HasScalar(variable))
                {
                    if (isdouble) this->AddDoubleScalar(variable, 0.0);
                    else this->AddIntegerScalar(variable, 0);
                }
                auto& v = *this->variableObjects[variable];
                if (isdouble != (v.GetType()==SurfaceVarType::Double))
                {
                    std::cout << "Setting incompatible types in SetScalarValues!" << std::endl;
                    abort();
                }
                for (size_t i = 0; i < this->numFaces; i++)
                {
                    v3d<double> c;
                    this->GetCentroid(i, &c[0], &c[1], &c[2]);
                    auto val = func(c[0], c[1], c[2]);
                    v.SetElem(i, val);
                }
            }
            
            void Rotate(v3d<double> axis, v3d<double> point, double theta)
            {
                auto rotateTrans = [&](double* x, double* y, double* z) -> void
                {
                    v3d<double> xyz(*x, *y, *z);
                    xyz -= point;
                    xyz.Rotate(axis, theta);
                    xyz += point;
                    *x = xyz[0];
                    *y = xyz[1];
                    *z = xyz[2];
                };
                
                this->Transform(rotateTrans);
            }
            
            AnalyticalGeometry& operator += (const AnalyticalGeometry& rhs)
            {
                double* newPoints = (double*)malloc(3*sizeof(double)*(this->numPoints + rhs.numPoints));
                int* newFaces = (int*)malloc(dimension*sizeof(int)*(this->numFaces + rhs.numFaces));
                
                for (int i = 0; i < dimension*this->numFaces; i++) newFaces[i] = this->faces[i];
                for (int i = 0; i < dimension*rhs.numFaces; i++) newFaces[i+dimension*this->numFaces] = this->numPoints+rhs.faces[i];
                
                for (int i = 0; i < 3*this->numPoints; i++) newPoints[i] = this->points[i];
                for (int i = 0; i < 3*rhs.numPoints; i++) newPoints[i+3*this->numPoints] = rhs.points[i];
                
                free(this->faces);
                free(this->points);
                
                this->faces = newFaces;
                this->points = newPoints;
                
                this->numFaces += rhs.numFaces;
                this->numPoints += rhs.numPoints;
                
                pidx = 3*numPoints;
                fidx = dimension*numFaces;
                
                return *this;
            }
            
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
            int numTimesAddedFace;
            std::map<std::string, double*> doubleScalars;
            std::map<std::string, int*> integerScalars;
            std::map<std::string, SurfaceVar*> variableObjects;
    };
}
#endif
