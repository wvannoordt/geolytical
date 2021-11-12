#pragma once

namespace geolytical
{
    /// @brief A simple 3D vector class
    /// \pre Note: use this class sparingly, as it is not optimized
    /// @author WVN
    template <typename numericType = double> struct v3d
    {
        /// @brief Zero constructor
        /// @author WVN
        v3d(void) {v[0] = 0; v[1] = 0; v[2] = 0;}
        
        /// @brief Constructor
        /// @param x X component
        /// @param y Y component
        /// @param z Z component
        /// @author WVN
        v3d(numericType x, numericType y, numericType z) {v[0] = x; v[1] = y; v[2] = z;}
        
        /// @brief Constructor
        /// @param x Pointer (of size > 3)
        /// @author WVN
        v3d(numericType* x) {v[0] = x[0]; v[1] = x[1]; v[2] = x[2];}
        
        /// @brief Constructor
        /// @param x value to initialize all elements with
        /// @author WVN
        template <typename othertype> v3d(othertype x) {v[0] = (othertype)x; v[1] = (othertype)x; v[2] = (othertype)x;}
        
        /// @brief Constructor, using initial and final values
        /// @param ini Pointer (of size > 3) to the coordinates of initial point
        /// @param ter Pointer (of size > 3) to the coordinates of terminal point
        /// @author WVN
        v3d(numericType* ini, numericType* ter) {v[0] = ter[0]-ini[0]; v[1] = ter[1]-ini[1]; v[2] = ter[2]-ini[2];}
        
        /// @brief Copy constructor
        /// @param w Vector to copy
        /// @author WVN
        v3d(const v3d& w) {v[0] = w.v[0]; v[1] = w.v[1]; v[2] = w.v[2];}
        
        /// @brief index operator
        /// @param i index
        /// @author WVN
        numericType & operator [] (int i) {return *(v+i);}
        
        ///@brief the data
        numericType v[3];
        
        ///@brief Normalize in-place by L2-norm
        void Normalize(void)
        {
            numericType norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            v[0]/=norm;
            v[1]/=norm;
            v[2]/=norm;
        }
        
        /// @brief Returns a string representing the vector
        std::string str()
        {
            std::string output = "[";
            output += std::to_string(v[0]) + ", ";
            output += std::to_string(v[1]) + ", ";
            output += std::to_string(v[2]) + "]";
            return output;
        }
        
        ///@brief Addition operator
        v3d operator + (const v3d& w) {return v3d(v[0]+w.v[0], v[1]+w.v[1], v[2]+w.v[2]);}
        
        v3d& operator += (const v3d& w) {v[0] += w.v[0]; v[1] += w.v[1]; v[2] += w.v[2]; return *this;}
        v3d& operator -= (const v3d& w) {v[0] -= w.v[0]; v[1] -= w.v[1]; v[2] -= w.v[2]; return *this;}
        
        
        ///@brief Difference operator
        v3d operator - (const v3d& w) {return v3d(v[0]-w.v[0], v[1]-w.v[1], v[2]-w.v[2]);}
        
        ///@brief Vector cross-product, a%b returns a "cross" b
        v3d operator % (const v3d& w) {return v3d(v[1]*w.v[2] - v[2]*w.v[1], v[2]*w.v[0] - v[0]*w.v[2], v[0]*w.v[1] - v[1]*w.v[0]);}
        
        ///@brief Scalar multiplication
        v3d operator * (const numericType a) {return v3d(a*v[0], a*v[1], a*v[2]);}
        
        v3d& operator *= (const numericType a) {v[0] *= a; v[1] *= a; v[2] *= a; return *this;}
        
        ///@brief Scalar division
        v3d operator / (const numericType a) {return v3d(v[0]/a, v[1]/a, v[2]/a);}
        
        v3d& Rotate(v3d& axis, double angle)
        {
            double cs = cos(angle);
            double ss = sin(angle);
            double ux = axis[0];
            double uy = axis[1];
            double uz = axis[2];
            double v1 = v[0];
            double v2 = v[1];
            double v3 = v[2];
            
            v[0] = (cs+ux*ux*(1.0-cs))*v1 + (ux*uy*(1.0-cs)-uz*ss)*v2 + (ux*uz*(1.0-cs)+uy*ss)*v3;
            v[1] = (uy*ux*(1.0-cs)+uz*ss)*v1 + (cs+uy*uy*(1.0-cs))*v2 + (uy*uz*(1.0-cs)-ux*ss)*v3;
            v[2] = (uz*ux*(1.0-cs)-uy*ss)*v1 + (uz*uy*(1.0-cs)+ux*ss)*v2 + (cs+uz*uz*(1.0-cs))*v3;
            return *this;
        }
        
        ///@brief L2 inner product
        numericType operator * (const v3d& w) {return v[0]*w.v[0] + v[1]*w.v[1] + v[2]*w.v[2];}
        
        v3d& RotateInZ(double angle)
        {
            double t1, t2;
            t1 = v[0];
            t2 = v[1];
            v[0] = t1*cos(angle) - t2*sin(angle);
            v[1] = t1*sin(angle) + t2*cos(angle);
            return *this;
        }
        
        ///@brief L2 norm
        numericType Norm(void) {return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);}
        
        v3d& operator = (const numericType rhs)
        {
            v[0] = rhs;
            v[1] = rhs;
            v[2] = rhs;
            return *this;
        }
    };
    
}