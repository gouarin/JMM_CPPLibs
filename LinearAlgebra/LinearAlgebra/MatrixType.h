//
//  LinearTransform.h
//  AmongConvex2
//
//  Created by Jean-Marie Mirebeau on 23/07/13.
//  Copyright (c) 2013 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef AmongConvex2_LinearTransform_h
#define AmongConvex2_LinearTransform_h

#include "VectorType.h"

namespace LinearAlgebra {
    
template<typename TComponent, size_t VRows, size_t VColumns>
struct Matrix :
vector_space< Matrix<TComponent,VRows,VColumns>, TComponent>
{
    typedef TComponent ComponentType;
    static const size_t Rows = VRows, Columns = VColumns;
    
    static bool IsInRange(int i, int j)  {return 0<=i && i<Rows && 0<=j && j<Columns;}
    ComponentType           & operator()(size_t i, size_t j)       {assert(IsInRange(int(i),int(j))); return data[LinearizedIndex(i,j)];}
    const ComponentType     & operator()(size_t i, size_t j) const {assert(IsInRange(int(i),int(j))); return data[LinearizedIndex(i,j)];}
    
    Matrix & operator+=(const Matrix & m){data+=m.data; return *this;}
    Matrix & operator-=(const Matrix & m){data-=m.data; return *this;}
    Matrix & operator*=(const ComponentType & a){data*=a;   return *this;}
    Matrix & operator/=(const ComponentType & a){data/=a;   return *this;}
    Matrix operator -() {Matrix m; m.data = -data; return m;}
    
    typedef Vector<ComponentType,VRows>    OutputVectorType;
    typedef Vector<ComponentType,VColumns>  InputVectorType;
    OutputVectorType operator * (const InputVectorType & u) const {
        OutputVectorType v;
        v.fill(ComponentType(0));
        for(size_t i=0; i<Rows; ++i){
            for(size_t j=0; j<Columns; ++j)
                v[i]+=this->operator()(i,j)*u[j];
        }
        return v;
    }
    
    template <size_t Columns2>
    Matrix<ComponentType, Rows, Columns2> operator * (const Matrix<ComponentType,Columns,Columns2> & m) const {
        Matrix<ComponentType, Rows, Columns2> p;
        p.fill(ComponentType(0));
        for(size_t i=0; i<Rows; ++i)
            for(size_t k=0; k<Columns2; ++k)
                for(size_t j=0; j<Columns; ++j)
                    p(i,k) += (this->operator()(i,j))*m(j,k);
        return p;
    }
    
    Matrix<ComponentType,Columns,Rows> Transpose() const {
        Matrix<ComponentType,Columns,Rows> m;
        for(size_t i=0; i<Rows; ++i)
            for(size_t j=0; j<Columns; ++j)
                m(j,i) = this->operator()(i,j);
        return m;
    }
    
    ComponentType Determinant() const;
    ComponentType Trace() const {
        static_assert(Rows==Columns,"Matrix must be square");
        ComponentType s(0);
        for(size_t i=0; i<Rows; ++i) s+=this->operator()(i,i);
        return s;
    }
    
    Matrix Inverse() const;
    ComponentType FrobeniusSquaredNorm() const {return data.SquaredNorm();}
    ComponentType FrobeniusNorm() const {return data.Norm();}
    
    void fill(const ComponentType & a){data.fill(a);}
    
    static Matrix Identity(){
        static_assert(Rows==Columns,"Matrix must be square");
        Matrix m;
        m.fill(ComponentType(0));
        for(size_t i=0; i<Rows; ++i) m(i,i)=1;
        return m;
    }
    
    // Do not forget double braces, like {{ {{1,2}} , {{3,4}} }}.
    static Matrix FromRows(const std::array<std::array<ComponentType, Columns>,Rows> & t){
        Matrix m;
        for(int i=0; i<Rows; ++i)
            for(int j=0; j<Columns; ++j)
                m(i,j) = t[i][j];
        return m;
    }
    
    static Matrix FromColumns(const std::array<std::array<ComponentType,Rows>,Columns> & t){
        Matrix m;
        for(int i=0; i<Rows; ++i)
            for(int j=0; j<Columns; ++j)
                m(i,j) = t[j][i];
        return m;
    }
    
    static Matrix Rotation(ComponentType theta){
        static_assert(Rows==2 && Columns==2,"Rotation matrices are two dimensional");
        static_assert( ! std::is_integral<ComponentType>::value,"Rotation matrices have real coefficients");
        const double c=cos(theta), s=sin(theta);
        return FromRows({{ {{c,-s}} , {{s,c}} }});
    }
    
    Vector<ComponentType,Rows*Columns> data;
protected:
    static size_t LinearizedIndex(size_t i, size_t j) {return i+j*Rows;}
};

template<typename TC, size_t VR, size_t VC>
typename Matrix<TC,VR,VC>::ComponentType
Matrix<TC,VR,VC>::Determinant() const
{
    static_assert(VR==VC,"Matrix must be square");
    static_assert(VR<=3,"Dimension <=3 only is supported");
    
    const Matrix & m = *this;
    switch (VR) {
        case 1: return m(0,0);
        case 2: return m(0,0)*m(1,1)-m(1,0)*m(0,1);
        case 3: {
            ComponentType det=0;
            for(int i=0; i<3; ++i) det+=m(i,0)*m((i+1)%3,1)*m((i+2)%3,2) - m(i,2)*m((i+1)%3,1)*m((i+2)%3,0);
            return det;
        }
        default: assert(false);
    }
}

template<typename TC, size_t VR, size_t VC>
Matrix<TC,VR,VC>
Matrix<TC,VR,VC>::Inverse() const
{
    static_assert(VR==VC,"Matrix must be square");
    static_assert(VR<=3,"Dimensions <=3 only are supported");
    Matrix m;
    const ComponentType d=Determinant();
    assert(d!=0); // exception ?
    switch (VR) {
        case 1: m(0,0)=1; break;
        case 2:
            m(0,0)= this->operator()(1,1);
            m(1,0)=-this->operator()(1,0);
            m(0,1)=-this->operator()(0,1);
            m(1,1)= this->operator()(0,0);
            break;
        case 3:
            for(int i=0; i<3; ++i)
                for(int j=0; j<3; ++j)
                    m(j,i)=
                    this->operator()((i+1)%3,(j+1)%3)*this->operator()((i+2)%3,(j+2)%3)-
                    this->operator()((i+1)%3,(j+2)%3)*this->operator()((i+2)%3,(j+1)%3);
            break;
        default: assert(false);
    }
    return m/d;
}
    
// Row based representation.
template<typename TC, size_t VR, size_t VC>
std::ostream & operator << (std::ostream & f, const Matrix<TC,VR,VC> & m)
{
    f<<"{";
    for(int i=0; i<VR; ++i){
        if(i>0) f<<",";
        f<<"{";
        for(int j=0; j<VC; ++j){
            if(j>0) f<<",";
            f<<m(i,j);
        }
        f<<"}";
    }
    f<<"}";
    return f;
}

}
#endif
