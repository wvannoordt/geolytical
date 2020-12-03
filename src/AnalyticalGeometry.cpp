#include "AnalyticalGeometry.h"
namespace geolytical
{
    AnalyticalGeometry::AnalyticalGeometry(void)
    {
        dealloc = false;
        dimension = -999;
        numPoints=-1;
        numFaces=-1;
        fidx=0;
        pidx=0;
        doComponentId=false;
        hasAnyScalars = false;
    }
    AnalyticalGeometry::~AnalyticalGeometry(void)
    {
        if (dealloc)
        {
            free(points);
            free(faces);
            free(compID);
            for (std::map<std::string, double*>::iterator it = doubleScalars.begin(); it!=doubleScalars.end(); it++)
            {
                free(it->second);
            }
            for (std::map<std::string, int*>::iterator it = integerScalars.begin(); it!=integerScalars.end(); it++)
            {
                free(it->second);
            }
            for (std::map<std::string, SurfaceVar*>::iterator it = variableObjects.begin(); it!=variableObjects.end(); it++)
            {
                delete it->second;
            }
        }
    }
    bool AnalyticalGeometry::HasScalar(std::string name)
    {
        return (doubleScalars.find(name)!=doubleScalars.end()) || (integerScalars.find(name)!=integerScalars.end());
    }
    void AnalyticalGeometry::Allocate(void)
    {
        dealloc = true;
        if (dimension<0) {std::cout << "dimension has not been set!" << std::endl; abort();}
        if (numPoints<0) {std::cout << "numPoints has not been set!" << std::endl; abort();}
        if (numFaces<0) {std::cout << "numFaces has not been set!" << std::endl; abort();}
        points = (double*)malloc(3*numPoints*sizeof(double));
        if (dimension==2) {faces = (int*)malloc(2*numFaces*sizeof(double));}
        else if (dimension==3) {faces = (int*)malloc(3*numFaces*sizeof(double));}
        else {std::cout << "dimension is neither 2 nor 3!" << std::endl; abort();}
        compID = (int*)malloc(numFaces*sizeof(int));
        for (int i = 0; i < numFaces; i++) compID[i] = 1;
    }
    void AnalyticalGeometry::AddPoint(double x, double y)
    {
        CHECKP;
        points[pidx] = x;
        pidx++;
        CHECKP;
        points[pidx] = y;
        pidx++;
    }
    void AnalyticalGeometry::AddPoint(double x, double y, double z)
    {
        CHECKP;
        points[pidx] = x;
        pidx++;
        CHECKP;
        points[pidx] = y;
        pidx++;
        CHECKP;
        points[pidx] = z;
        pidx++;
    }
    void AnalyticalGeometry::AddFace(int x, int y)
    {
        CHECKF;
        faces[fidx] = x;
        fidx++;
        CHECKF;
        faces[fidx] = y;
        fidx++;
    }
    void AnalyticalGeometry::AddFace(int x, int y, int z)
    {
        CHECKF;
        faces[fidx] = x;
        fidx++;
        CHECKF;
        faces[fidx] = y;
        fidx++;
        CHECKF;
        faces[fidx] = z;
        fidx++;
    }
    void AnalyticalGeometry::AddFaceSwitch(int x, int y, int z)
    {
        CHECKF;
        faces[fidx] = y;
        fidx++;
        CHECKF;
        faces[fidx] = x;
        fidx++;
        CHECKF;
        faces[fidx] = z;
        fidx++;
    }
    void AnalyticalGeometry::AddFace(int x, int y, int z, bool switchValues)
    {
        if (switchValues)
        {
            AddFaceSwitch(x, y, z);
        }
        else
        {
            AddFace(x, y, z);
        }
    }
    void AnalyticalGeometry::ErrIfOBPoint(void)
    {
        if (pidx>=3*numPoints) {std::cout << "point index out of bounds!" << std::endl; abort();}
    }
    void AnalyticalGeometry::ErrIfOBFace(void)
    {
        if (fidx>=dimension*numFaces) {std::cout << "face index out of bounds!" << std::endl; abort();}
    }
    void AnalyticalGeometry::ErrIfNotAllPoint(void)
    {
        if (pidx != 3*numPoints) {std::cout << "output vtk before all points filled (" << pidx << "/" << 3*numPoints << ")" << std::endl; abort();}
    }
    void AnalyticalGeometry::ErrIfNotAllFace(void)
    {
        if (fidx != dimension*numFaces) {std::cout << "output vtk before all faces filled (" << fidx << "/" << dimension*numFaces << ")" << std::endl; abort();}
    }
    void AnalyticalGeometry::GetCentroid(int faceIdx, double* xout, double* yout, double* zout)
    {
        int p1 = faces[faceIdx*dimension+0];
        int p2 = faces[faceIdx*dimension+1];
        int p3 = faces[faceIdx*dimension+dimension-1];
        double x1 = points[3*p1+0];
        double x2 = points[3*p2+0];
        double x3 = points[3*p3+0]*(dimension-2);
        double y1 = points[3*p1+1];
        double y2 = points[3*p2+1];
        double y3 = points[3*p3+1]*(dimension-2);
        double z1 = points[3*p1+2];
        double z2 = points[3*p2+2];
        double z3 = points[3*p3+2]*(dimension-2);
        *xout = (1.0/dimension) * (x1+x2+x3);
        *yout = (1.0/dimension) * (y1+y2+y3);
        *zout = (1.0/dimension) * (z1+z2+z3);
    }
    
    void AnalyticalGeometry::RemoveScalar(std::string name)
    {
        if (HasScalar(name))
        {
            SurfaceVar* obj = variableObjects[name];
            if (obj->GetType() == SurfaceVarType::Integer)
            {
                free(integerScalars[name]);
                integerScalars.erase(name);
                std::cout << name << std::endl;
            }
            else
            {
                free(doubleScalars[name]);
                doubleScalars.erase(name);
            }
            delete variableObjects[name];
            variableObjects.erase(name);
        }
    }
    
    SurfaceVar& AnalyticalGeometry::AddDoubleScalar(std::string name, double value)
    {
        if (HasScalar(name)) {std::cout << "Attempt to add duplicate scalar \"" << name << "\"!" << std::endl; abort();}
        double* array = (double*)malloc(numFaces*sizeof(double));
        for (int i = 0; i < numFaces; i++) array[i] = value;
        doubleScalars.insert({name, array});
        hasAnyScalars = true;
        variableObjects.insert({name, new SurfaceVar(name, this, array, numFaces, SurfaceVarType::Double)});
        return *(variableObjects[name]);
    }

    SurfaceVar& AnalyticalGeometry::AddDoubleScalar(std::string name, double(*func)(int))
    {
        if (HasScalar(name)) {std::cout << "Attempt to add duplicate scalar \"" << name << "\"!" << std::endl; abort();}
        double* array = (double*)malloc(numFaces*sizeof(double));
        for (int i = 0; i < numFaces; i++) array[i] = func(i);
        doubleScalars.insert({name, array});
        hasAnyScalars = true;
        variableObjects.insert({name, new SurfaceVar(name, this, array, numFaces, SurfaceVarType::Double)});
        return *(variableObjects[name]);
    }

    SurfaceVar& AnalyticalGeometry::AddDoubleScalar(std::string name, double(*func)(double, double, double))
    {
        double x, y, z;
        if (HasScalar(name)) {std::cout << "Attempt to add duplicate scalar \"" << name << "\"!" << std::endl; abort();}
        double* array = (double*)malloc(numFaces*sizeof(double));
        for (int i = 0; i < numFaces; i++)
        {
            GetCentroid(i, &x, &y, &z);
            array[i] = func(x, y, z);
        }
        doubleScalars.insert({name, array});
        hasAnyScalars = true;
        variableObjects.insert({name, new SurfaceVar(name, this, array, numFaces, SurfaceVarType::Double)});
        return *(variableObjects[name]);
    }

    SurfaceVar& AnalyticalGeometry::AddIntegerScalar(std::string name, int value)
    {
        if (HasScalar(name)) {std::cout << "Attempt to add duplicate scalar \"" << name << "\"!" << std::endl; abort();}
        int* array = (int*)malloc(numFaces*sizeof(int));
        for (int i = 0; i < numFaces; i++) array[i] = value;
        integerScalars.insert({name, array});
        hasAnyScalars = true;
        variableObjects.insert({name, new SurfaceVar(name, this, array, numFaces, SurfaceVarType::Integer)});
        return *(variableObjects[name]);
    }

    SurfaceVar& AnalyticalGeometry::AddIntegerScalar(std::string name, int (*func)(int))
    {
        if (HasScalar(name)) {std::cout << "Attempt to add duplicate scalar \"" << name << "\"!" << std::endl; abort();}
        int* array = (int*)malloc(numFaces*sizeof(int));
        for (int i = 0; i < numFaces; i++) array[i] = func(i);
        integerScalars.insert({name, array});
        hasAnyScalars = true;
        variableObjects.insert({name, new SurfaceVar(name, this, array, numFaces, SurfaceVarType::Integer)});
        return *(variableObjects[name]);
    }
    
    SurfaceVar& AnalyticalGeometry::AddIntegerScalar(std::string name, int (*func)(double, double, double))
    {
        double x, y, z;
        if (HasScalar(name)) {std::cout << "Attempt to add duplicate scalar \"" << name << "\"!" << std::endl; abort();}
        int* array = (int*)malloc(numFaces*sizeof(int));
        for (int i = 0; i < numFaces; i++)
        {
            GetCentroid(i, &x, &y, &z);
            array[i] = func(x, y, z);
        }
        integerScalars.insert({name, array});
        hasAnyScalars = true;
        variableObjects.insert({name, new SurfaceVar(name, this, array, numFaces, SurfaceVarType::Integer)});
        return *(variableObjects[name]);
    }

    void AnalyticalGeometry::OutputToVtk(std::string filename)
    {
        FULLP;
        FULLF;
        if (doComponentId) AddIntegerScalar("Components", 1);
        std::ofstream myfile;
        myfile.open(filename.c_str());
        myfile << "# vtk DataFile Version 3.0" << std::endl;
        myfile << "vtk output" << std::endl;
        myfile << "ASCII" << std::endl;
        myfile << "DATASET POLYDATA" << std::endl;
        myfile << "POINTS " << numPoints << " double" << std::endl;
        for (int i = 0; i < numPoints; i++)
        {
            myfile << points[3*i] << " " << points[3*i+1] << " " << points[3*i+2] << std::endl;
        }
        myfile << ((this->dimension==3)?("POLYGONS "):("LINES ")) << numFaces << " " << (1+dimension)*numFaces << std::endl;
        for (int i = 0; i < numFaces; i++)
        {
            myfile << dimension;
            for (int d = 0; d < dimension; d++)
            {
                myfile << " " << faces[i*dimension+d];
            }
            myfile << std::endl;
        }
        if (hasAnyScalars)
        {
            myfile << "CELL_DATA " << numFaces << std::endl;
            for (std::map<std::string, double*>::iterator it = doubleScalars.begin(); it!=doubleScalars.end(); it++)
            {
                myfile << "SCALARS " << it->first << " double" << std::endl;
                myfile << "LOOKUP_TABLE default" << std::endl;
                for (int i = 0; i < numFaces; i++) myfile << (it->second)[i] << std::endl;
            }
            for (std::map<std::string, int*>::iterator it = integerScalars.begin(); it!=integerScalars.end(); it++)
            {
                myfile << "SCALARS " << it->first << " int" << std::endl;
                myfile << "LOOKUP_TABLE default" << std::endl;
                for (int i = 0; i < numFaces; i++) myfile << (it->second)[i] << std::endl;
            }
        }
        myfile.close();
    }
    void AnalyticalGeometry::Deform(Transformation3D deformer)
    {
        for (int i = 0; i < numPoints; i++)
        {
            deformer(points+3*i, points+3*i+1, points+3*i+2);
        }
    }
    void AnalyticalGeometry::Deform(Transformation2D deformer)
    {
        for (int i = 0; i < numPoints; i++)
        {
            deformer(points+3*i, points+3*i+1);
        }
    }
    void AnalyticalGeometry::OutputPointsAsCSV(std::string filename)
    {
        FULLP;
        std::ofstream myfile;
        myfile.open(filename.c_str());
        for (int i = 0; i < numPoints; i++)
        {
            myfile << points[3*i] << ", " << points[3*i+1] << ", " << points[3*i+2] << std::endl;
        }
        myfile.close();
    }
}