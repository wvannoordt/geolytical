#ifndef GEO_VTK_TRIANGU_H
#define GEO_VTK_TRIANGU_H
#include "AnalyticalGeometry.h"
namespace geolytical
{
    class VtkTriangulation : public AnalyticalGeometry
    {
        public:
            VtkTriangulation(std::string vtkFilename_in);
        private:
            std::string vtkFilename;
            std::string header;
            std::string metaString;
            bool ReadFile(std::string file, std::string& errMessage);
    };
}

#endif