//Description of analytical function describing the shape of a symmetrical 4-digit NACA-00xx airfoil
double getAirfoilY_NACA00xx(int& xx, double x)
{
  double t = (double) xx/100.0;
  double y = 5*t*(0.2969*sqrt(x) - 0.1260*x - 0.3516*pow(x,2) + 0.2843*pow(x,3) - 0.1036*pow(x,4));
  
  return y;
}

//Mirrors the point distribution for a NACA-00xx symmetrical airfoil across the chord and
//rearranges data in a structure suitable for the VTK writer
//Assumes the data for the upper-half of the airfoil is given
std::vector<double> mirrorAndRearrangeNACA00xxData(std::vector<double>& xy_top)
{
  int size_xy = xy_top.size();
  int nxpts   = size_xy/2;

  int totalPts = 4*nxpts - 4;
  std::vector<double> xy(totalPts,0.0);  
  for (int i = 0; i < nxpts; i++)
    {
      if (i == 0)
        {
          xy[0] = 1.0;  //x
          xy[1] = 0.0;  //y
        }
      else if (i == (nxpts-1))
        {
          xy[2*i  ] = 0.0;  //x
          xy[2*i+1] = 0.0;  //y
        }
      else
        {
          xy[2*(i-1)+2]          =  xy_top[size_xy - 2 - 2*i    ];   //  x
          xy[2*(i-1)+3]          = -xy_top[size_xy - 2 - 2*i + 1];   // -y
          xy[totalPts-2*(i-1)-2] =  xy_top[size_xy - 2 - 2*i    ];   //  x
          xy[totalPts-2*(i-1)-1] =  xy_top[size_xy - 2 - 2*i + 1];   // +y
        }
    }
  return xy;
}

//Creates a uniformly spaced out distribution of  grid points across the whole chord length
//based on the number of points
std::vector<double> uniformlyDistributedNACA00xx(int& xx, int& nxpts, double& dx)
{
  int totalPts = 4*nxpts-4;
  std::vector<double> xy(totalPts,0.0);
  for (int i = 0; i < nxpts; i++)
    {
      double x = (double) i/(nxpts-1);     
      if (i == 0)
        {
          xy[0] = 1.0;  //x
          xy[1] = 0.0;  //y
        }
      else if (i == (nxpts-1))
        {
          xy[2*i  ] = 0.0;  //x
          xy[2*i+1] = 0.0;  //y
        }
      else
        {
          xy[2*(i-1)+2]          =  1.0 - x;                    // x
          xy[2*(i-1)+3]          = -getAirfoilY_NACA00xx(xx,x); //-y
          xy[totalPts-2*(i-1)-2] =  1.0 - x;                    // x
          xy[totalPts-2*(i-1)-1] =  getAirfoilY_NACA00xx(xx,x); //+y
        }
    }
  
  return xy;
}

//Creates a point distribution for a NACA-00xx airfoil between two arbitraty x values
std::vector<double> NACA00xxPointDistribution(int& xx, int& npts, double& xstart, double& xend)
{
  double dx = (xend - xstart)/(npts-1);
  std::vector<double> xy;
  for (int i = 0; i < npts; i++)
    {
      double x = (double) xstart + i*dx;
      double y = getAirfoilY_NACA00xx(xx,x);
      xy.push_back(x);
      xy.push_back(y);
    }

  return xy;
}

//Concatenates point distributions from two separate regions given that the last point of the first
//region and the first point of the second region is the same
std::vector<double> concatenatePointDistributions(std::vector<double>& xy1, std::vector<double>& xy2)
{
   int n_xy1 = xy1.size();
   int n_xy2 = xy2.size();

  //Check if (x,y) for the last point in xy1 is the same as (x,y) for the first point in xy2
   double epsilon = 1.0e-8;
   double x1 = xy1[n_xy1-2];   double y1 = xy1[n_xy1-1];
   double x2 = xy2[0];         double y2 = xy2[1];

//   std::cout << x1  << " " << x2 << std::endl;
//   std::cout << y1  << " " << y2 << std::endl;
//   std::cout << x1 - x2 << " " << y1 - y2 << std::endl;
   
   if ((std::abs(x1-x2) > epsilon) || (std::abs(y1-y2) > epsilon))
     {
       std::cout << "[E] Continuity of arrays broken." << std::endl;
       std::cout << "[E] Last element of first vector doesn't match the first element of the second vector." << std::endl;
       std::cout << "[E] Code exiting." << std::endl;
       abort();      
     }

   int n_xy = n_xy1+n_xy2-2;
   std::vector<double> xy(n_xy,0.0);
   for (int i = 0; i < (n_xy1 - 2); i++)
     {
       xy[i] = xy1[i];
     }   
   
   for (int i = (n_xy1 - 2); i < n_xy; i++)
     {
       xy[i] = xy2[i - n_xy1 + 2];
     }

   return xy;
}
