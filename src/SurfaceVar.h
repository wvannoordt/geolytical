#ifndef GEO_DURF_VAR_H
#define GEO_DURF_VAR_H
#include <iostream>
#include <string>
namespace geolytical
{
    class AnalyticalGeometry;
    namespace SurfaceVarType
    {
        enum SurfaceVarType
        {
            Double,
            Integer
        };
    }
    class SurfaceVar
    {
        public:
            SurfaceVar(std::string name_in, AnalyticalGeometry* domain_in, void* array_in, size_t num_in, SurfaceVarType::SurfaceVarType type_in);
            ~SurfaceVar(void);
            SurfaceVarType::SurfaceVarType GetType(void) {return type;}
            void SetElem(size_t i, double val)
            {
                double* p = (double*)array;
                p[i] = val;
            }
            void SetElem(size_t i, int    val)
            {
                int* p = (int*)array;
                p[i] = val;
            }
            void* Arr(void) {return array;}
            size_t GetSize(void) {return num;}
        private:
            void Build(std::string name_in, AnalyticalGeometry* domain_in, void* array_in, size_t num_in, SurfaceVarType::SurfaceVarType type_in);
            std::string name;
            AnalyticalGeometry* domain;
            void* array;
            size_t num;
            SurfaceVarType::SurfaceVarType type;
    };
}

#endif