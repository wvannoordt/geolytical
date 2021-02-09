#include "GeoTypes.h"
#include "Cylinder.h"
#include "DeformationTypes.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <cmath>
#include <unistd.h>
namespace geolytical
{
    Cylinder::Cylinder(int nx_in_cyl, int nr_in_rad, int layersonFace_in, bbox bounds_in)
    {
        dimension = 3;
        bounds = bounds_in;
        nx = nx_in_cyl+1;
        nr = nr_in_rad;
        layersonFace = layersonFace_in;
        halvingFactor = 0;
        frequencyOfHalving = 100000;
        numberOfRadialPointsAtLevel = new int[layersonFace];
        offx = bounds.xmin;
        offy = (bounds.ymax + bounds.ymin) * 0.5;
        offz = (bounds.zmax + bounds.zmin) * 0.5;
        scalex = bounds.xmax - bounds.xmin;
        scaley = (bounds.ymax - bounds.ymin) * 0.5;
        scalez = (bounds.zmax - bounds.zmin) * 0.5;
        for (int i = 0; i < 8*sizeof(int); i++)
        {
            if ((layersonFace & (1<<i)) == 0)
            {
                halvingFactor++;
            }
            else break;
        }
        pointIdx = 0;
        numberOfRadialPointsAtLevel[layersonFace-1] = nr;
        int factor2level = 1;
        double minAcceptable = 0.1;
        for (int i = 0; i < layersonFace-1; i++)
        {
            int level = layersonFace - 2 - i;
            numberOfRadialPointsAtLevel[layersonFace - 2 - i] = numberOfRadialPointsAtLevel[layersonFace - 2 - (i - 1)];
            int ni = numberOfRadialPointsAtLevel[layersonFace - 2 - i];
            double dsi = ((double)level/(double)layersonFace) / ni;
            if (dsi < minAcceptable/nr)
            {
                if (numberOfRadialPointsAtLevel[layersonFace - 2 - i]%2!=0)
                {
                    std::cout << "geolytical::Cylinder: nr_in_rad not divisible enough by 2!" << std::endl;
                    abort();
                }
                ni = ni/2;
                numberOfRadialPointsAtLevel[layersonFace - 2 - i] = ni;
            }
        }
        int numPointAtInner = nr_in_rad / layersonFace;
        isIrregularLayer = new bool[layersonFace-1];
        for (int i = 0; i < layersonFace-1; i++)
        {
            isIrregularLayer[i] = numberOfRadialPointsAtLevel[i]!=numberOfRadialPointsAtLevel[i+1];
        }
        CountPoints();
        CountFaces();
        Allocate();
        CreatePoints();
        WriteFaceTriangles();
        WriteCylTriangles();
    }
    
    int Cylinder::GetCylinderRingStart(int ilayer)
    {
        if (ilayer == 0) return pointsOnFace - nr;
        if (ilayer == nx-1) return 2*pointsOnFace - nr;
        return 2*pointsOnFace + (ilayer-1)*nr;
    }
    
    void Cylinder::CountFaces(void)
    {
        numFaces = (nx-1)*nr;
        numFaces += numberOfRadialPointsAtLevel[0];//core face
        for (int i = 1; i < layersonFace; i++)
        {
            if (!isIrregularLayer[i-1])
            {
                numFaces += 2*numberOfRadialPointsAtLevel[i];
            }
            else
            {
                numFaces += 3*numberOfRadialPointsAtLevel[i]/2;
            }
        }
        numFaces*= 2;
    }
    
    void Cylinder::CreatePoints(void)
    {
        CreateFacePoints();
        CreateCylPoints();
    }
    
    void Cylinder::CreateFacePoints(void)
    {
        double dr = 1.0/layersonFace;
        for (int zface = 0; zface <= 1; zface++)
        {
            WritePoint(0, 0, (double)zface);
            for (int ilayer = 0; ilayer < layersonFace; ilayer++)
            {
                double dtheta = 6.28318530718 / numberOfRadialPointsAtLevel[ilayer];
                for(int j = 0; j < numberOfRadialPointsAtLevel[ilayer]; j++)
                {
                    double theta = j * dtheta;
                    double r = (ilayer+1)*dr;
                    double x = r*cos(theta);
                    double y = r*sin(theta);
                    double z = (double)zface;
                    WritePoint(x, y, z);
                }
            }
        }
    }
    
    void Cylinder::WritePoint(double x, double y, double z)
    {
        //switch z and x
        AddPoint(offx + scalex*z, offy + scaley*y, offz + scalez*x);
    }
    
    void Cylinder::WriteFace(bool reverseOrder, int p1, int p2, int p3)
    {
        AddFace(p1, p2, p3, reverseOrder);
    }
    
    void Cylinder::CreateCylPoints(void)
    {
        double dtheta = 6.28318530718 / nr;
        double dz = 1.0/(nx-1);
        for (int ix = 1; ix < nx-1; ix++)
        {
            for (int ir = 0; ir < nr; ir++)
            {
                double theta = ir*dtheta;
                double z = ix*dz;
                WritePoint(cos(theta), sin(theta), z);
            }
        }
    }
    
    void Cylinder::CountPoints(void)
    {
        int numPointsOnCyl = (nx-2) * nr; // exclude boundaries
        int numPointsOnFace = 1; //axis point
        for (int i = 0; i < layersonFace; i++)
        {
            numPointsOnFace += numberOfRadialPointsAtLevel[i];
        }
        numPoints = numPointsOnCyl + 2*numPointsOnFace;
        pointsOnFace = numPointsOnFace;
    }
    
    Cylinder::~Cylinder(void)
    {
        if (dealloc)
        {
            dealloc = false;
            delete [] numberOfRadialPointsAtLevel;
            delete [] isIrregularLayer;
            free(points);
        }
    }
    
    int Cylinder::GetLayerStart(int ilayer)
    {
        int offset = 0;
        if (ilayer == 0) return 1 + offset;
        else
        {
            int output = 1;
            for (int i = 0; i < ilayer; i++)
            {
                output += numberOfRadialPointsAtLevel[i];
            }
            return output + offset;
        }
    }
    
    void Cylinder::WriteFaceTriangles(void)
    {
        for (int z = 0; z <=1; z++)
        {
            int totalOffset = 0;
            if (z==1) totalOffset = pointsOnFace;
            for (int i = 0; i < numberOfRadialPointsAtLevel[0]; i++)
            {
                WriteFace(z==1, totalOffset + 0, totalOffset + 1+i, totalOffset + 1 + ((1+i)%numberOfRadialPointsAtLevel[0]));
            }
            for (int layer = 0; layer < layersonFace-1; layer++)
            {
                int numPointHere = numberOfRadialPointsAtLevel[layer];
                if (isIrregularLayer[layer])
                {
                    int start1 = GetLayerStart(layer);
                    int start2 = GetLayerStart(layer+1);
                    for (int i = 0; i < numberOfRadialPointsAtLevel[layer]; i++)
                    {
                        WriteFace(z!=1, totalOffset + start1 + i, totalOffset + start1 + (i + 1)%numPointHere, totalOffset + start2 + 2*i+1);
                        WriteFace(z!=1, totalOffset + start1 + i, totalOffset + start2 + 2*i+1, totalOffset + start2 + 2*i);
                        WriteFace(z!=1, totalOffset + start2 + 2*i+1, totalOffset + start1 + (i+1)%numPointHere, totalOffset + start2 + 2*((i+1)%numPointHere));
                    }
                }
                else
                {
                    int start1 = GetLayerStart(layer);
                    int start2 = GetLayerStart(layer+1);
                    for (int i = 0; i < numberOfRadialPointsAtLevel[layer]; i++)
                    {
                        WriteFace(z!=1, totalOffset + start1 + i, totalOffset + start1 + (i + 1)%numPointHere, totalOffset + start2 + i);
                        WriteFace(z!=1, totalOffset + start1 + (i + 1)%numPointHere, totalOffset + start2 + (i + 1)%numPointHere, totalOffset + start2 + i);
                    }
                }
            }
        }
    }
    
    void Cylinder::WriteCylTriangles(void)
    {
        for (int i = 0; i < nx-1; i++)
        {
            int start1 = GetCylinderRingStart(i);
            int start2 = GetCylinderRingStart(i+1);
            for (int j = 0; j < nr; j++)
            {
                WriteFace(true, start1 + j, start1 + (j + 1)%nr, start2 + j);
                WriteFace(true, start1 + (j + 1)%nr, start2 + (j + 1)%nr, start2 + j);
            }
        }
    }
}
