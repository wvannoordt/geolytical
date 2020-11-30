#include "SurfaceVar.h"
namespace geolytical
{
    SurfaceVar::SurfaceVar(std::string name_in, AnalyticalGeometry* domain_in, void* array_in, size_t num_in, SurfaceVarType::SurfaceVarType type_in)
    {
        Build(name_in, domain_in, array_in, num_in, type_in);
    }
    
    SurfaceVar::~SurfaceVar(void)
    {
        
    }
    
    void SurfaceVar::Build(std::string name_in, AnalyticalGeometry* domain_in, void* array_in, size_t num_in, SurfaceVarType::SurfaceVarType type_in)
    {
        name_in = name;
        domain_in = domain;
        array_in = array;
        num_in = num;
        type_in = type;
    }
}