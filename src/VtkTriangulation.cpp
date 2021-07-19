#include "VtkTriangulation.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>
#include <string>
namespace geolytical
{
    VtkTriangulation::VtkTriangulation(std::string vtkFilename_in)
    {
        vtkFilename = vtkFilename_in;
        std::string errMessage;
        if (!ReadFile(vtkFilename, errMessage))
        {
            std::cout << errMessage << std::endl;
            abort();
        }
    }
    
    bool VtkTriangulation::ReadFile(std::string file, std::string& errMessage)
    {
        dimension = 3;
        std::ifstream myfile;
        myfile.open(file);
        std::string line;
        std::string junk;
        std::string token;
        errMessage = "error parsing header material";
        auto next = [&](std::string& ll) -> decltype(std::getline(myfile, line)) {return std::getline(myfile, ll);};
        if (!next(header)) return false;
        if (!next(metaString)) return false;
        if (!next(line)) return false;
        if (line.find("ASCII")==std::string::npos)
        {
            errMessage = "Only ascii files are supported!";
            return false;
        }
        if (!next(line)) return false;
        token = "DATASET POLYDATA";
        errMessage = "Expecting token \"" + token +"\" line: " + line;
        if (line.find(token)==std::string::npos)
        {
            return false;
        }
        
        if (!next(line))
        {
            errMessage = "Cannot find info about \"POINTS\"";
            return false;
        }
        errMessage = "Error parsing line \"" + line + "\"";
        std::istringstream buffer(line);
        buffer >> junk;
        buffer >> numPoints;
        buffer >> junk;
        if (buffer.fail()) return false;
        AllocatePoints();
        
        double x, y, z;
        for (int i = 0; i < numPoints; i++)
        {
            if (!next(line))
            {
                errMessage = "Not enough points specified";
                return false;
            }
            std::istringstream ptbuf(line);
            ptbuf >> x;
            ptbuf >> y;
            ptbuf >> z;
            if (buffer.fail())
            {
                errMessage = "Bad point on the following line: " + line;
                return false;
            }
            AddPoint(x, y, z);
        }
        
        if (!next(line))
        {
            errMessage = "Cannot find info about \"POLYGONS\"";
            return false;
        }
        errMessage = "Error parsing line \"" + line + "\"";
        std::istringstream faceBuf(line);
        int dummy;
        faceBuf >> junk;
        faceBuf >> numFaces;
        faceBuf >> dummy;
        if (faceBuf.fail()) return false;
        AllocateFaces();
        int ct, n1, n2, n3;
        for (int i = 0; i < numFaces; i++)
        {
            if (!next(line))
            {
                errMessage = "Not enough faces specified";
                return false;
            }
            errMessage = "Error parsing line \"" + line + "\"";
            std::istringstream faceCtBuf(line);
            faceCtBuf >> ct;
            faceCtBuf >> n1;
            faceCtBuf >> n2;
            faceCtBuf >> n3;
            if (ct != 3)
            {
                errMessage = "Only triangular faces are supported!";
                return false;
            }
            if (faceCtBuf.fail()) return false;
            AddFace(n1, n2, n3);
        }
        
        while (line.find("CELL_DATA")==std::string::npos && next(line)){}
        while(next(line))
        {
            errMessage = "Error parsing scalar data header";
            std::istringstream namebuf(line);
            std::string name, varType;
            namebuf >> junk;
            namebuf >> name;
            namebuf >> varType;
            
            if (!next(line)) return false;
            token = "LOOKUP_TABLE default";
            errMessage = "Expecting token \"" + token + "\" line: " + line;
            if (line.find(token)==std::string::npos)
            {
                return false;
            }
            if (varType == "double")
            {
                auto& var = this->AddDoubleScalar(name, 0.0);
                errMessage = "not enough elements for scalar \"" + name + "\"";
                double* arr = (double*)var.Arr();
                double n;
                for (int i = 0; i < numFaces; i++)
                {
                    if (!next(line)) return false;
                    std::istringstream db(line);
                    db >> n;
                    if (db.fail())
                    {
                        errMessage = "Cannot parse line \"" + line + "\" in scalar \"" + name + "\"";
                        return false;
                    }
                    arr[i] = n;
                }
            }
            else if (varType == "int")
            {
                auto& var = this->AddIntegerScalar(name, 0);
                int* arr = (int*)var.Arr();
                errMessage = "not enough elements for scalar \"" + name + "\"";
                int n;
                for (int i = 0; i < numFaces; i++)
                {
                    if (!next(line)) return false;
                    std::istringstream db(line);
                    db >> n;
                    if (db.fail())
                    {
                        errMessage = "Cannot parse line \"" + line + "\" in scalar \"" + name + "\"";
                        return false;
                    }
                    arr[i] = n;
                }
            }
            else
            {
                errMessage = "Only supporting scalar data of type \"int\" and \"double\"";
                return false;
            }
        }
        errMessage = "No Error";
        return true;
    }
}