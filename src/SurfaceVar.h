#ifndef GEO_DURF_VAR_H
#define GEO_DURF_VAR_H

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