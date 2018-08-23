// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef LiftedFastMarching_SymmetricMatrixType_h
#define LiftedFastMarching_SymmetricMatrixType_h

#include "MatrixType.h"

namespace LinearAlgebra {
    
template<typename TComponent, size_t VDimension>
struct SymmetricMatrix :
vector_space< SymmetricMatrix<TComponent, VDimension>, TComponent>
{
    typedef TComponent ComponentType;
    static const size_t Dimension = VDimension;
    typedef Vector<ComponentType, Dimension> VectorType;
    typedef Point<ComponentType, Dimension> PointType;
    typedef Matrix<ComponentType, Dimension, Dimension> MatrixType;
    
    // Coefficient access
    static bool IsInRange(int i, int j){
        return 0<=i && i<Dimension && 0<=j && j<Dimension;}
    static int LinearizedIndex(int i, int j);
    ComponentType & operator()(int i, int j){return data[LinearizedIndex(i,j)];}
    const ComponentType & operator()(int i, int j) const {return data[LinearizedIndex(i,j)];}

    // Linear operations
    SymmetricMatrix & operator+=(const SymmetricMatrix & m){data+=m.data; return *this;}
    SymmetricMatrix & operator-=(const SymmetricMatrix & m){data-=m.data; return *this;}
    SymmetricMatrix & operator*=(const ComponentType & a){data*=a; return *this;}
    SymmetricMatrix & operator/=(const ComponentType & a){data/=a; return *this;}
    SymmetricMatrix operator -(){SymmetricMatrix m; m.data=-data; return m;}
    
    // Geometry
    template<typename T> ComponentType
    ScalarProduct(const Vector<T, Dimension> &, const Vector<T, Dimension> &) const;
    
    template<typename T> ComponentType
    SquaredNorm(const Vector<T,Dimension> & u) const {return ScalarProduct(u, u);}
    
    template<typename T> ComponentType
    Norm(const Vector<T, Dimension> & u) const {return sqrt(SquaredNorm(u));}
    
    VectorType Gradient(const VectorType & u) const {return operator*(u)/Norm(u);}
    
    template<typename T> bool
    IsAcute(const Vector<T, Dimension> & u, const Vector<T, Dimension> & v) const {
        return ScalarProduct(u, v)>=0;}

    template<typename T> VectorType
    operator*(const Vector<T, Dimension> &) const;
    
    ComponentType Trace() const;
    ComponentType Determinant() const;
    SymmetricMatrix Inverse() const;
    VectorType CGSolve(VectorType) const;
    
    // Constructors
    static SymmetricMatrix Zero();
    static SymmetricMatrix Identity();
    static SymmetricMatrix Diagonal(const VectorType &);
    template<typename T> static SymmetricMatrix RankOneTensor(const Vector<T, Dimension> & u);
    static SymmetricMatrix RandomPositive();
    static SymmetricMatrix FromUpperTriangle(const MatrixType & mat);
    SymmetricMatrix operator=(ComponentType a); // Dimension 1 only
    
    template<typename T, size_t D> SymmetricMatrix<ComponentType,D>
    Gram(const std::array<Vector<T,Dimension>, D> & a) const;
    template<size_t D> SymmetricMatrix<ComponentType, D> // a^T.m.a
    Gram(const Matrix<ComponentType,Dimension,D> & a) const;
    template<size_t D> SymmetricMatrix<ComponentType, D> // a.m.a^T
    GramT(const Matrix<ComponentType,D,Dimension> & a) const;
    template<typename T, size_t D>
    static SymmetricMatrix EuclideanGram(const std::array<Vector<T,D>, Dimension> & a);

    // Conversion
    explicit operator MatrixType() const;
    template<typename T> static SymmetricMatrix CastCoordinates(const SymmetricMatrix<T, Dimension> & m0);
    
    static const size_t InternalDimension = (Dimension*(Dimension+1))/2;
    Vector<ComponentType,InternalDimension> data;
    template<typename ...T,typename dummy = typename std::enable_if<sizeof...(T)==InternalDimension>::type >
    constexpr SymmetricMatrix(T... t):data(t...){};

    SymmetricMatrix(){};
protected:
    const ComponentType & coef(int i, int j) const {return this->operator()(i,j);}
};


    // Printing
    template<typename TC, size_t VD>
    std::ostream & operator << (std::ostream & f, const SymmetricMatrix<TC,VD> & m)
    {
        f<<"{";
        for(int i=0; i<VD; ++i){
            if(i>0) f<<",";
            f<<"{";
            for(int j=0; j<VD; ++j){
                if(j>0) f<<",";
                f<<m(i,j);
            }
            f<<"}";
        }
        f<<"}";
        return f;
    }

#include "Implementation/SymmetricMatrixType.hxx"
    
} // end namespace LinearAlgebra

#endif
