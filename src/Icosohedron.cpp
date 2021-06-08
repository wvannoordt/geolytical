#include "Icosohedron.h"
#include <cmath>
#include "v3d.h"
#include <vector>
#define PI 3.141592653589793238462643383279
namespace geolytical
{
    Icosohedron::Icosohedron(int nxEquator_in, bbox bounds_in)
    {
        nFaceOnEdge = this->GetFacesOnEdge(nxEquator_in);
        numTrisOnFace = this->GetNumTrisOnFace(nFaceOnEdge);
        bounds = bounds_in;
        dimension = 3;
        CountFaces();
        CountPoints();
        Allocate();
        CreatePoints();
        CreateFaces();
        bbox boundsTemp = bounds;
        bounds.xmin = -1;
        bounds.xmax =  1;
        bounds.ymin = -1;
        bounds.ymax =  1;
        bounds.zmin = -1;
        bounds.zmax =  1;
        this->RemapBoundingBox(boundsTemp);
    }

    int Icosohedron::GetFacesOnEdge(int nxEquator)
    {
        return nxEquator;
    }

    int Icosohedron::GetNumTrisOnFace(int nFacesOnEdge)
    {
        return ((nFacesOnEdge)*(nFacesOnEdge-1)/2) + (nFacesOnEdge*(nFacesOnEdge+1)/2);
    }

    void Icosohedron::CountFaces(void)
    {
        numFaces = 20*numTrisOnFace;
    }

    void Icosohedron::CountPoints(void)
    {
        int numBaseVertices = 12;
        int numPointsOnEdge = nFaceOnEdge + 1;
        int numEdges = 30;
        int numPointsOnBaseVertex = 1;
        int numPointsOnBaseEdges = numPointsOnEdge*numEdges;
        int numBaseFaces = 20;
        int numPointsOnFace = (nFaceOnEdge+1)*(nFaceOnEdge+2)/2;
        numPoints = numBaseFaces*numPointsOnFace - numEdges*numPointsOnEdge + numBaseVertices*numPointsOnBaseVertex;
    }

    void Icosohedron::CreatePoints(void)
    {
        double phi = 0.5 + 0.5*sqrt(5.0);
        double phibar = 0.5 - 0.5*sqrt(5.0);
        AddPoint(0.0, 0.0, -1.0);
        double edgeLength = 1.0/(sin(2.0*(PI)/5.0));
        double dx = edgeLength/nFaceOnEdge;
        const double a = sqrt(3.0)/2.0;
        const double b = 1.0;
        const double c = 0.5*sqrt(5.0 + 2.0*sqrt(5.0));
        const double angleOfEdgeFromVertexNormal = acos((c*c + b*b - a*a) / (2.0*b*c));
        const double angleOfFaceFromVertexNormal = acos((c*c + a*a - b*b) / (2.0*a*c));
        const double dihedralAngle = acos(-sqrt(5.0)/3.0);
        const double pentAngle = 3.0*PI/5.0;
        const double pentAngleComp = 2.0*PI/5.0;
        v3d<> x(0.0, 0.0, -1);
        v3d<> dxVecEdge(dx*cos(angleOfEdgeFromVertexNormal), 0.0, dx*sin(angleOfEdgeFromVertexNormal));
        v3d<> dxFace(-dx, 0.0, 0.0);
        dxFace.RotateInZ(0.5*pentAngle);
        for (int i = 0; i < nFaceOnEdge; i++)
        {
            x += dxVecEdge;
            for (int edg = 0; edg < 5; edg++)
            {
                for (int j = 0; j <= i; j++)
                {
                    x += dxFace;
                    AddPoint(x[0], x[1], x[2]);
                }
                dxFace.RotateInZ(-pentAngleComp);
            }
        }
        dxVecEdge[0] = 0.0;
        dxVecEdge[1] = 0.0;
        dxVecEdge[2] = 1.0;
        dxFace[0] = 0.0;
        dxFace[1] = -1.0;
        dxFace[2] = 0.0;
        dxFace.RotateInZ(-0.5*pentAngleComp+0.5*PI);
        dxVecEdge.Rotate(dxFace, PI/6.0);
        dxFace.RotateInZ(-0.5*PI);
        dxVecEdge.Rotate(dxFace, angleOfFaceFromVertexNormal+0.5*PI-dihedralAngle);
        dxFace*=dx;
        dxVecEdge *= dx;
        for (int i = 0; i < nFaceOnEdge-1; i++)
        {
            x += dxVecEdge;
            for (int edg = 0; edg < 5; edg++)
            {
                for (int j = 0; j < nFaceOnEdge-1-i; j++)
                {
                    x += dxFace;
                    AddPoint(x[0], x[1], x[2]);
                }
                dxFace.RotateInZ(-0.5*pentAngleComp);
                for (int j = 0; j <= i; j++)
                {
                    x += dxFace;
                    AddPoint(x[0], x[1], x[2]);
                }
                dxFace.RotateInZ(-0.5*pentAngleComp);
            }
        }
        x += dxVecEdge;
        dxFace.RotateInZ(-0.5*pentAngleComp);
        dxVecEdge[0] = 0.0;
        dxVecEdge[1] = 0.0;
        dxVecEdge[2] = dx;
        v3d<> axis(0.0, 1.0, 0.0);
        axis.RotateInZ(0.5*pentAngle-0.5*PI);
        dxVecEdge.Rotate(axis, angleOfEdgeFromVertexNormal-0.5*PI);
        for (int i = 0; i < nFaceOnEdge; i++)
        {
            if (i!=0) x += dxVecEdge;
            for (int edg = 0; edg < 5; edg++)
            {
                for (int j = 0; j < nFaceOnEdge-i; j++)
                {
                    x += dxFace;
                    AddPoint(x[0], x[1], x[2]);
                }
                dxFace.RotateInZ(-pentAngleComp);
            }
        }
        AddPoint(0.0, 0.0, 1.0);
    }

    void Icosohedron::CreateFaces(void)
    {
        std::vector<int> layerCount;
        layerCount.push_back(1);
        for (int i = 1; i <= nFaceOnEdge; i++)
        {
            layerCount.push_back(5*i);
        }
        int midLayer = layerCount[layerCount.size()-1];
        for (int i = 0; i < nFaceOnEdge-1; i++)
        {
            layerCount.push_back(midLayer);
        }
        for (int i = nFaceOnEdge; i >= 1; i--)
        {
            layerCount.push_back(5*i);
        }
        layerCount.push_back(1);
        std::vector<int> layerBegin;
        layerBegin.reserve(layerCount.size());
        layerBegin.push_back(0);
        for (int i = 0; i < layerCount.size()-1; i++) layerBegin.push_back(layerBegin[i] + layerCount[i]);
        int end = numPoints-1;
        for (int i = 0; i < 5; i++)
        {
            AddFace(0, (i%5)+1, ((i+1)%5)+1);
        }
        for (int i = 1; i <= nFaceOnEdge-1; i++)
        {
            int rootLo = layerBegin[i];
            int nlo = layerCount[i];
            int rootHi = layerBegin[i+1];
            int nhi = layerCount[i+1];
            auto numbrHi = [=](int k) -> int{return rootHi+(k+nhi-1)%nhi;};
            auto numbrLo = [=](int k) -> int{return rootLo+(k+nlo-1)%nlo;};
            for (int edg = 0; edg < 5; edg++)
            {
                AddFace(numbrLo(edg*i), numbrHi(edg*(i+1)), numbrHi(1+edg*(i+1)));
                for (int j = 0; j < i; j++)
                {
                    AddFace(numbrLo(edg*i+j+1), numbrLo(edg*i+j), numbrHi(edg*(i+1)+j+1));
                    AddFace(numbrLo(edg*i+j+1), numbrHi(edg*(i+1)+j+1), numbrHi(edg*(i+1)+j+2));
                }
            }
        }
        
        for (int i = 1; i <= nFaceOnEdge; i++)
        {
            int rootLo = layerBegin[i+nFaceOnEdge-1];
            int nlo = layerCount[i+nFaceOnEdge-1];
            int rootHi = layerBegin[i+nFaceOnEdge];
            int nhi = layerCount[i+nFaceOnEdge];
            auto numbrHi = [=](int k) -> int{return rootHi+(k+nhi-1)%nhi;};
            auto numbrLo = [=](int k) -> int{return rootLo+(k+nlo-1)%nlo;};
            for (int j = 0; j < nlo; j++)
            {
                AddFace(numbrLo(j), numbrHi(j), numbrLo(j+1));
                AddFace(numbrLo(j+1), numbrHi(j), numbrHi(j+1));
            }
        }
        
        for (int i = 0; i < nFaceOnEdge-1; i++)
        {
            int rootLo = layerBegin[i+2*nFaceOnEdge];
            int nlo = layerCount[i+2*nFaceOnEdge];
            int rootHi = layerBegin[i+2*nFaceOnEdge+1];
            int nhi = layerCount[i+2*nFaceOnEdge+1];
            int ieff = nFaceOnEdge-1 - i;
            auto numbrHi = [=](int k) -> int{return rootHi+(k+nhi-1)%nhi;};
            auto numbrLo = [=](int k) -> int{return rootLo+(k+nlo-1)%nlo;};
            for (int edg = 0; edg < 5; edg++)
            {
                AddFace(numbrLo(edg*(ieff+1)), numbrHi(edg*ieff), numbrLo(edg*(ieff+1)+1));
                for (int j = 0; j < ieff; j++)
                {
                    AddFace(numbrLo(edg*(ieff+1)+j+1), numbrHi(edg*ieff + j), numbrHi(edg*ieff + j + 1));
                    AddFace(numbrLo(edg*(ieff+1)+j+1), numbrHi(edg*ieff + j + 1), numbrLo(edg*(ieff+1)+j+2));
                }
            }
        }
        for (int i = 0; i < 5; i++)
        {
            AddFace(end, end-(i%5)-1, end-((i+1)%5)-1);
        }
    }

    Icosohedron::~Icosohedron(void)
    {
        
    }
}