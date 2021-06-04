#include "GeoTypes.h"
#include "FlatPlate.h"
#include "DeformationTypes.h"
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
namespace geolytical
{
    FlatPlate::FlatPlate(int nx_in, int nz_in, bbox bounds_in)
    {
        dimension = 3;
        bounds = bounds_in;
        deleteLevels = false;
        nx = nx_in + 1;
        nz = nz_in + 1;
        numLevels = 0;
        enableAlignedCenters = false;
        CountPoints();
        CountFaces();
        Allocate();
        CreatePoints();
        CreateFaces();
    }
    
    void FlatPlate::CountPoints(void)
    {
        int numLevelsx, numLevelsz;
        numLevelsx = 0;
        numLevelsz = 0;
        int numx = nx;
        int numz = nz;
        while (numx>1)
        {
            numLevelsx++;
            numx = 1+((numx-1)/2);
        }
        while (numz>1)
        {
            numLevelsz++;
            numz = 1+((numz-1)/2);
        }
        numLevels = GEO_MAX(numLevelsx, numLevelsz);
        numx = nx;
        numz = nz;
        levelx = new int[numLevels];
        levelz = new int[numLevels];
        numPointsInLayer = new int[numLevels];
        upperLineArray = new int[GEO_MAX(nx,nz)];
        lowerLineArray = new int[GEO_MAX(nx,nz)];
        minPointIdArray = new int[GEO_MAX(nx,nz)];
        levelPtStart = new int[numLevels];
        deleteLevels = true;
        for (int i = numLevels-1; i >= 0; i--)
        {
            levelx[i] = GEO_MAX(2,numx);
            levelz[i] = GEO_MAX(2,numz);
            numx = 1+((numx-1)/2);
            numz = 1+((numz-1)/2);
        }
        int totalPoints = 0;
        levelPtStart[0] = 0;
        for (int i = 0; i < numLevels; i++)
        {
            int layerTotalPoints = (2*levelx[i] + 2*levelz[i] - 4);
            if (i == (numLevels-1))
            {
                layerTotalPoints = levelx[i]*levelz[i];
            }
            numPointsInLayer[i] = layerTotalPoints;
            totalPoints += layerTotalPoints;
            if (i>0)
            {
                levelPtStart[i] = levelPtStart[i-1] + numPointsInLayer[i-1];
            }
        }
        numPoints = totalPoints;
    }
    
    void FlatPlate::CountFaces(void)
    {
        int facesOnBottom = 2;
        int facesOnTop = 2*(nx-1)*(nz-1);
        int facesZ = 0;
        int facesX = 0;
        for (int i = 1; i < numLevels; i++)
        {
            facesZ += (levelz[i] + levelz[i-1] - 2);
            facesX += (levelx[i] + levelx[i-1] - 2);
        }
        numFaces = facesOnBottom + facesOnTop + 2*(facesZ + facesX);
    }
    
    void FlatPlate::GetLineArray(bool isX, bool upperFace, int level)
    {
        int startIdxLow = levelPtStart[level];
        int startIdxHigh = levelPtStart[level+1];
        int offset = 0;
        if (isX)
        {
            offset = (upperFace?levelx[level]:0);
            nLower = levelx[level];
            for (int i = 0; i < nLower; i++) lowerLineArray[i] = startIdxLow + offset + i;
            offset = (upperFace?levelx[level+1]:0);
            nUpper = levelx[level+1];
            for (int i = 0; i < nUpper; i++) upperLineArray[i] = startIdxHigh + offset + i;
        }
        else
        {
            nLower = levelz[level];
            nUpper = levelz[level+1];
            offset = 2*levelx[level+1] - 1 + (upperFace?(levelz[level+1]-2):(0));
            for (int i = 1; i < nUpper-1; i++)
            {
                upperLineArray[i] = startIdxHigh + offset + i;
            }
            offset = 2*levelx[level] - 1 + (upperFace?(levelz[level]-2):(0));
            for (int i = 1; i < nLower-1; i++)
            {
                lowerLineArray[i] = startIdxLow + offset + i;
            }
            lowerLineArray[0] = startIdxLow  + (upperFace?(levelx[level]-1):0);
            upperLineArray[0] = startIdxHigh + (upperFace?(levelx[level+1]-1):0);
            lowerLineArray[nLower-1] = startIdxLow  + (upperFace?(2*levelx[level]-1):levelx[level]);
            upperLineArray[nUpper-1] = startIdxHigh + (upperFace?(2*levelx[level+1]-1):levelx[level+1]);
        }
        
        if (level == numLevels-2)
        {
            int offsetTop = 0;
            int pitch = isX?1:nx;
            if (upperFace)
            {
                if (isX)
                {
                    offsetTop = (nz-1)*nx;
                }
                else
                {
                    offsetTop = nx-1;
                }
            }
            int startIdx = startIdxHigh + offsetTop;
            for (int i = 0; i < nUpper; i++)
            {
                upperLineArray[i] = startIdx+i*pitch;
            }
        }
    }
    
    void FlatPlate::SetEnableAlignedCenters(bool val)
    {
        enableAlignedCenters = val;
        ResetFaceCounter();
        CountFaces();
        CreateFaces();
    }
    
    void FlatPlate::CreateFaces(void)
    {
        int n = 0;
        AddFace(0, 1, 2);
        AddFace(2, 1, 3);
        
        for (int upperLower = 0; upperLower < 2; upperLower++)
        {
            std::string ulstring = ((upperLower==1)?"upper":"lower");
            for (int Xthenz = 0; Xthenz < 2; Xthenz++)
            {
                std::string dstring = ((Xthenz==1)?"z":"x");
                for (int level = 0; level < numLevels-1; level++)
                {
                    GetLineArray((Xthenz==0), (upperLower==1), level);
                    if ((nUpper == 2) && (nLower == 2))
                    {
                        AddFace(upperLineArray[1], upperLineArray[0], lowerLineArray[0], (upperLower==1)==(Xthenz==1));
                        AddFace(upperLineArray[1], lowerLineArray[0], lowerLineArray[1], (upperLower==1)==(Xthenz==1));
                    }
                    else
                    {
                        for (int i = 0; i < nUpper-1; i++)
                        {
                            int lowerTarget = (int)round(nLower*((double)i/nUpper));
                            if (i==0) lowerTarget = 0;
                            if (i==nUpper-2) lowerTarget = nLower-1;                    
                            AddFace(upperLineArray[i+1], upperLineArray[i], lowerLineArray[lowerTarget], (upperLower==1)==(Xthenz==1));
                            minPointIdArray[lowerTarget] = i+1;
                        }
                        for (int i = 0; i < nLower-1; i++)
                        {
                            AddFace(lowerLineArray[i], lowerLineArray[i+1], upperLineArray[minPointIdArray[i]], (upperLower==1)==(Xthenz==1));
                        }
                    }
                }
            }
        }
        if (enableAlignedCenters)
        {
            for (int ix = 0; ix <  nx-1; ix++)
            {
                int p1x, p2x, p3x;
                p1x = levelPtStart[numLevels-1] + (0)*(nx) + (ix);
                p2x = levelPtStart[numLevels-1] + (0)*(nx) + (ix+1);
                p3x = levelPtStart[numLevels-1] + (1)*(nx) + (ix+1);
                AddFace(p1x, p3x, p2x);
                for (int iz = 0; iz < nz-2; iz++)
                {
                    int p1,p2,p3,p4;
                    p1 = levelPtStart[numLevels-1] + (iz+0)*(nx) + (ix+0);
                    p2 = levelPtStart[numLevels-1] + (iz+1)*(nx) + (ix+0);
                    p3 = levelPtStart[numLevels-1] + (iz+1)*(nx) + (ix+1);
                    p4 = levelPtStart[numLevels-1] + (iz+2)*(nx) + (ix+1);
                    AddFace(p4, p1, p2);
                    AddFace(p3, p1, p4);
                }
                p1x = levelPtStart[numLevels-1] + (nz-2)*(nx) + (ix);
                p2x = levelPtStart[numLevels-1] + (nz-1)*(nx) + (ix);
                p3x = levelPtStart[numLevels-1] + (nz-1)*(nx) + (ix+1);
                AddFace(p1x, p2x, p3x);
            }
        }
        else
        {
            for (int iz = 0; iz < nz-1; iz++)
            {
                for (int ix = 0; ix < nx-1; ix++)
                {
                    int p1,p2,p3,p4;
                    p1 = levelPtStart[numLevels-1] + (iz+0)*(nx) + (ix+0);
                    p2 = levelPtStart[numLevels-1] + (iz+0)*(nx) + (ix+1);
                    p3 = levelPtStart[numLevels-1] + (iz+1)*(nx) + (ix+0);
                    p4 = levelPtStart[numLevels-1] + (iz+1)*(nx) + (ix+1);
                    AddFace(p2, p1, p3);
                    AddFace(p4, p2, p3);
                }
            }
        }
    }
    
    void FlatPlate::CreatePoints(void)
    {
        int numWritten = 0;
        double dy = (bounds.ymax - bounds.ymin)/(numLevels-1);
        for (int i = 0; i < numLevels-1; i++)
        {
            int nxHere = levelx[i];
            int nzHere = levelz[i];
            double dx = (bounds.xmax - bounds.xmin) / (nxHere - 1);
            double dz = (bounds.zmax - bounds.zmin) / (nzHere - 1);
            for (int faceIdx = 0; faceIdx < 2; faceIdx++)
            {
                for (int ix = 0; ix < nxHere; ix++)
                {
                    AddPoint(bounds.xmin + ix*dx, bounds.ymin+i*dy, faceIdx*bounds.zmax + (1-faceIdx)*bounds.zmin);
                    numWritten++;
                }
            }
            for (int faceIdx = 0; faceIdx < 2; faceIdx++)
            {
                for (int iz = 1; iz < nzHere-1; iz++)
                {
                    AddPoint(faceIdx*bounds.xmax + (1-faceIdx)*bounds.xmin, bounds.ymin+i*dy, bounds.zmin + iz*dz);
                    numWritten++;
                }
            }
        }
        for (int iz = 0; iz < nz; iz++)
        {
            double dz = (bounds.zmax - bounds.zmin)/(nz-1);
            for (int ix = 0; ix < nx; ix++)
            {
                double dx = (bounds.xmax - bounds.xmin)/(nx-1);
                AddPoint(bounds.xmin + ix*dx, bounds.ymax, bounds.zmin + iz*dz);
            }
        }
    }
    
    FlatPlate::~FlatPlate(void)
    {
        if (deleteLevels)
        {
            deleteLevels = false;
            delete [] levelx;
            delete [] levelz;
            delete [] numPointsInLayer;
            delete [] upperLineArray;
            delete [] lowerLineArray;
            delete [] minPointIdArray;
            delete [] levelPtStart;
        }
    }
}
