// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

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
    
    InputVectorType Solve(OutputVectorType) const;
    
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
    template<size_t d, typename Dummy> struct Det;
};

// Row based printing.
template<typename TC, size_t VR, size_t VC>
std::ostream & operator << (std::ostream & f, const Matrix<TC,VR,VC> & m) {
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
    
// -------- Determinant --------
template<typename TC, size_t VR, size_t VC> template<typename Dummy>
struct Matrix<TC,VR,VC>::Det<0, Dummy> {
    typedef Matrix<TC,VR,VC> M;
    M::ComponentType operator()(const M & m){return TC(1);}
};

template<typename TC, size_t VR, size_t VC> template<typename Dummy>
struct Matrix<TC,VR,VC>::Det<1, Dummy> {
    typedef Matrix<TC,VR,VC> M;
    M::ComponentType operator()(const M & m){return m(0,0);}
};

template<typename TC, size_t VR, size_t VC> template<typename Dummy>
struct Matrix<TC,VR,VC>::Det<2, Dummy> {
    typedef Matrix<TC,VR,VC> M;
    M::ComponentType operator()(const M & m){return m(0,0)*m(1,1)-m(1,0)*m(0,1);}
};
    
template<typename TC, size_t VR, size_t VC> template<typename Dummy>
struct Matrix<TC,VR,VC>::Det<3, Dummy> {
    typedef Matrix<TC,VR,VC> M;
    M::ComponentType operator()(const M & m){
        M::ComponentType det=0;
        for(int i=0; i<3; ++i) det+=m(i,0)*m((i+1)%3,1)*m((i+2)%3,2) - m(i,2)*m((i+1)%3,1)*m((i+2)%3,0);
        return det;
    }
};
    
template<typename TC, size_t VR, size_t VC> template<size_t n, typename Dummy>
struct Matrix<TC,VR,VC>::Det {
    typedef Matrix<TC,VR,VC> M;
    M::ComponentType operator()(const M & m){
        // Get largest coefficient, in absolute value
        assert(!std::numeric_limits<TC>::is_integer);
        const size_t d=VR;
        int iMax,jMax;
        M::ComponentType vMax = TC(0);
        for(int i=0; i<d; ++i)
            for(int j=0; j<d; ++j){
                const ComponentType v = std::abs(m(i,j));
                if(v>vMax){
                    iMax=i; jMax=j; vMax=v;
                }
            }
        if(vMax==TC(0)) return TC(0);
        vMax=m(iMax,jMax);
        
        // Pivot
        Matrix<TC, d-1, d-1> a; a.fill(TC(0));
        for(int i=0, ia=0; i<d; ++i){
            if(i==iMax) continue;
            const ComponentType delta = m(i,jMax)/vMax;
            for(int j=0, ja=0; j<d; ++j){
                if(j==jMax) continue;
                a(ia,ja)=m(i,j)-delta*m(iMax,j);
                ++ja;
            }
            ++ia;
        }
        return a.Determinant() *vMax *((iMax+jMax)%2==0 ? 1 : -1);
    }
};

template<typename TC, size_t VR, size_t VC>
typename Matrix<TC,VR,VC>::ComponentType
Matrix<TC,VR,VC>::Determinant() const
{
    static_assert(VR==VC,"Matrix must be square");
    return Det<VR,void>()(*this);
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
    
template<typename TC, size_t VR, size_t VC> auto
Matrix<TC,VR,VC>::Solve(OutputVectorType b) const -> InputVectorType {
    // A basic Gauss pivot
    static_assert(VR==VC,"Matrix must be square");
    const size_t n=VR;
    std::array<int,n> i2j, j2i; i2j.fill(-1); j2i.fill(-1);
    Matrix m = *this;
    for(int j=0; j<n; ++j){
        // Get largest coefficient in column
        ComponentType cMax=0;
        int iMax=0;
        for(int i=0; i<n; ++i){
            if(i2j[i]>=0) continue;
            const ComponentType c = m(i,j);
            if(std::abs(c)>std::abs(cMax)){
                cMax=c; iMax=i;}
        }
        i2j[iMax]=j;
        j2i[j]=iMax;
        assert(cMax!=0); // Matrix is not invertible
        
        // Remove line from other lines, while performing likewise on b
        for(int i=0; i<n; ++i){
            if(i2j[i]>=0) continue;
            const ComponentType r = m(i,j)/cMax;
            for(int k=j+1; k<n; ++k)
                m(i,k)-=m(iMax,k)*r;
            b[i]-=b[iMax]*r;
        }
    }
    // Solve remaining triangular system
    InputVectorType a;
    for(int j=n-1; j>=0; --j){
        const int i=j2i[j];
        ComponentType & r = a[j];
        r=b[i];
        for(int k=j+1; k<n; ++k){
            r-=a[k]*m(i,k);}
        r/=m(i,j);
    }
    return a;
}

    


}
#endif
