#define PI 3.14159265359
#include <iostream>
#include <vector>
#include <cmath>
#include "geolytical.h"
#include "naca.h"

int main()
{
  //Chord length (also acts as scaling factor if chord length is not unity)
  double chordLength = 1.0;

  //Description for a symmetric, 4-digit NACA airfoil
  //NACA-00xx
  int xx  = 12;
  
  //Number of regions with grid spacings other than the default grid spacing
  //nRefinementRegions = 0 gives a uniformly distributed VTK across the chord length
  int nRefinementRegions = 3;
 
  std::vector<double>  xRefinement(nRefinementRegions,0.1);

  int nRegions = nRefinementRegions+1;
  std::vector<int>  nPtRegions(nRegions,1);
  nPtRegions[0]  = 201;
  if (nRegions > 1)
    {
      nPtRegions[1]  = 101;
      nPtRegions[2]  = 51;
      nPtRegions[3]  = 101;
      xRefinement[0] = 0.05;
      xRefinement[1] = 0.20;
      xRefinement[2] = 0.75;
    }
  
  // Default grid-spacing for the airfoil assuming no refinement regions exist, i.e., uniform grid-spacing from x = 0 to 1
  int nxloc = nPtRegions[0];
  double defaultDx = 1.0/nxloc;

  //nxloc = Total number of points along the chord (includes all points across all regions, 0 & 1 included)
  std::vector<double> xloc(2*nxloc,0.0), yloc(2*nxloc,0.0);
  
  if (nRefinementRegions == 0)
    {
      std::vector<double> xy_data;
      xy_data = uniformlyDistributedNACA00xx(xx, nxloc, defaultDx);
      geolytical::ExtrudedCurve2D foil3D(10, 0.0, 0.2, xy_data.data(), xy_data.size()/2);
      foil3D.OutputToVtk("output3D.vtk");
      return 0;      
    }
  else
    {      
      double xstart = 0.0;
      double xend   = 0.0;
      int npts      = 0;
      std::vector<double> xy_0, xy_n, xy_i, xy_data, xy_data_final;
      for (int i = 0; i < nRegions; i++)
        {
          if (i == 0)
            {
              xstart   = 0.0;
              xend     = xRefinement[0];
              npts     = nPtRegions[0];
              xy_0     = NACA00xxPointDistribution(xx, npts, xstart, xend);
            }
          else if (i == (nRegions - 1))
            {
              xstart  = xRefinement.back();
              xend    = 1.0;
              npts    = nPtRegions.back();
              xy_n    = NACA00xxPointDistribution(xx, npts, xstart, xend);
              if (nRegions == 2)
                {
                  xy_data = concatenatePointDistributions(xy_0,xy_n);
                }
              else
                {
                  xy_data = concatenatePointDistributions(xy_data,xy_n);
                }
            }
          else
            {
              xstart   = xRefinement[i-1];
              xend     = xRefinement[i  ];
              npts     = nPtRegions[i];
              xy_i     = NACA00xxPointDistribution(xx, npts, xstart, xend);
              if (i == 1)
                {
                  xy_data  = concatenatePointDistributions(xy_0,xy_i);
                }
              else
                {
                  xy_data  = concatenatePointDistributions(xy_data,xy_i);
                }
            }
        } 
      xy_data_final = mirrorAndRearrangeNACA00xxData(xy_data);
      geolytical::ExtrudedCurve2D foil3D(21, 0.0, 0.2, xy_data_final.data(), xy_data_final.size()/2);
      foil3D.AddIntegerScalar("Components",[](double x, double y, double z)
                              {
                                if (x < 0.075) return 2;
                                if (x > 0.075 && x < 0.85 &&y > 0.0) return 3;
                                if (x >= 0.85) return 4;                                  
                                return 1;                                
                              });
      foil3D.OutputToVtk("output3D_2.vtk");      
      return 0;            
    }
}
