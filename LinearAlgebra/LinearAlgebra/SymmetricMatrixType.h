//
//  SymmetricMatrixType.h
//  LiftedFastMarching
//
//  Created by Jean-Marie Mirebeau on 09/10/2014.
//  Copyright (c) 2014 Jean-Marie Mirebeau. All rights reserved.
//

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
    
    static bool IsInRange(int i, int j){
        return 0<=i && i<Dimension && 0<=j && j<Dimension;}
    ComponentType & operator()(int i, int j){
        return data[LinearizedIndex(i,j)];}
    const ComponentType & operator()(int i, int j) const {
        return data[LinearizedIndex(i,j)];}

    SymmetricMatrix & operator+=(const SymmetricMatrix & m){data+=m.data; return *this;}
    SymmetricMatrix & operator-=(const SymmetricMatrix & m){data-=m.data; return *this;}
    SymmetricMatrix & operator*=(const ComponentType & a){data*=a; return *this;}
    SymmetricMatrix & operator/=(const ComponentType & a){data/=a; return *this;}
    SymmetricMatrix operator -(){SymmetricMatrix m; m.data=-data; return m;}
    
    template<typename T>
    ComponentType
    ScalarProduct(const Vector<T, Dimension> & u, const Vector<T, Dimension> & v)
    const {
        ComponentType sum=0;
        for(int i=0; i<Dimension; ++i){
            ComponentType sumi=0;
            for(int j=0; j<Dimension; ++j)
                sumi+=coef(i,j)*ComponentType(v[j]);
            sum+=ComponentType(u[i])*sumi;
        }
        return sum;
    }
    
    template<typename T>
    ComponentType SquaredNorm(const Vector<T,Dimension> & u) const {
        return ScalarProduct(u, u);
    }
    
    template<typename T>
    ComponentType Norm(const Vector<T, Dimension> & u) const {
        return sqrt(SquaredNorm(u));
    }
    
    template<typename T>
    VectorType operator*(const Vector<T, Dimension> & u) const {
        VectorType v;
        for(int i=0; i<Dimension; ++i){
            v[i]=0;
            for(int j=0; j<Dimension; ++j)
                v[i]+=coef(i, j)*ComponentType(u[j]);
        }
        return v;
    }
    
    ComponentType Determinant() const;
    SymmetricMatrix Inverse() const;
    
    static SymmetricMatrix Zero(){
        SymmetricMatrix m;
        std::fill(m.data.begin(), m.data.end(), ComponentType(0));
        return m;
    }
    
    static SymmetricMatrix Identity(){
        SymmetricMatrix m;
        for(int i=0; i<Dimension; ++i)
            for(int j=0; j<=i; ++j)
                m(i,j)=(i==j);
        return m;
    }
    
    template<typename T>
    static SymmetricMatrix RankOneTensor(const Vector<T, Dimension> & u){
        SymmetricMatrix m;
        for(int i=0; i<Dimension; ++i)
            for(int j=0; j<=i; ++j)
                m(i,j)=u[i]*u[j];
        return m;
    }
    
    static SymmetricMatrix RandomPositive(){
        MatrixType mat;
        for(int i=0; i<Dimension; ++i)
            for(int j=0; j<Dimension; ++j)
                mat(i,j) = -1+2*ComponentType(std::rand())/RAND_MAX;
        return FromUpperTriangle(mat.Transpose()*mat);
/*        MatrixType symMat = mat.Transpose()*mat;
        SymmetricMatrix m;
        for(int i=0; i<Dimension; ++i)
            for(int j=0; j<=i; ++j)
                m(i,j) = symMat(i,j);
        return m;*/
    }
    
    static SymmetricMatrix FromUpperTriangle(const MatrixType & mat){
        SymmetricMatrix m;
        for(int i=0; i<Dimension; ++i)
            for(int j=i; j<Dimension; ++j)
                m(i,j)=mat(i,j);
        return m;
    }
    
    SymmetricMatrix operator=(ComponentType a){
        static_assert(Dimension==1,"Assignment from scalar to matrix in dimension 1 only");
        data[0]=a;
        return *this;
    }
    
    static int LinearizedIndex(int i, int j){
        assert(IsInRange(i, j));
        if(i<j) std::swap(i, j);
        return (i*(i+1))/2+j;
    }
    
    ComponentType Trace() const {
        ComponentType result=0;
        for(int i=0; i<Dimension; ++i)
            result+=coef(i,i);
        return result;
    }
    
    explicit operator MatrixType(){
        MatrixType result;
        for(int i=0; i<Dimension; ++i)
            for(int j=0; j<Dimension; ++j)
                result(i,j)=coef(i, j);
        return result;
    }
    
    template<typename T, int D>
    SymmetricMatrix<ComponentType,D> Gram(const std::array<Vector<T,Dimension>, D> & a) const {
        SymmetricMatrix<ComponentType,D> m;
        for(int i=0; i<D; ++i)
            for(int j=0; j<=i; ++j)
                m(i,j) = ScalarProduct(a[i],a[j]);
        return m;
    }
    
    template<typename T, int D>
    static SymmetricMatrix EuclideanGram(const std::array<Vector<T,D>, Dimension> & a){
        SymmetricMatrix m;
        for(int i=0; i<Dimension; ++i)
            for(int j=0; j<=i; ++j)
                m(i,j) = a[i].ScalarProduct(a[j]);
        return m;
    }

    template<typename T>
    bool IsAcute(const Vector<T, Dimension> & u, const Vector<T, Dimension> & v) const {
        return ScalarProduct(u, v)>=0;}
    
protected:
    static const size_t InternalDimension = (Dimension*(Dimension+1))/2;
    Vector<ComponentType,InternalDimension> data;
    const ComponentType & coef(int i, int j) const {return this->operator()(i,j);}
};


    template<typename TC, size_t VD>
    TC
    SymmetricMatrix<TC, VD>::Determinant() const {
        static_assert(VD>0 && VD<=3,"Unsupported matrix size");
        switch (VD) {
            case 0: return 1;
            case 1: return coef(0,0);
            case 2: return coef(0,0)*coef(1,1)-coef(1,0)*coef(0,1);
            case 3: {
                ComponentType result=0;
                for(int i=0; i<3; ++i){
                    const int j=(i+1)%3, k=(i+2)%3;
                    result = result
                    +coef(i,0)*coef(j,1)*coef(k,2)
                    -coef(k,0)*coef(j,1)*coef(i,2);
                }
                return result;
            }
            default: throw "Unsupported matrix size";
        }
    }
    
    template<typename TC, size_t VD>
    SymmetricMatrix<TC, VD>
    SymmetricMatrix<TC, VD>::Inverse() const {
        ComponentType det = Determinant();
        SymmetricMatrix  m;
        switch (VD) {
            case 0: return SymmetricMatrix();
            case 1: m(0,0) = ComponentType(1)/det; return m;
            case 2: m(0,0) = coef(1,1); m(1,1) = coef(0,0); m(0,1)=-coef(0,1);
                m/=det; return m;
            case 3:
                for(int i=0; i<3; ++i)
                    for(int j=0; j<=i; ++j){
                        int i1 = (i+1)%3, i2=(i+2)%3, j1=(j+1)%3, j2=(j+2)%3;
                        m(i,j) = coef(i1,j1)*coef(i2,j2)-coef(i1,j2)*coef(i2,j1);
                    }
                m/=det; return m;
            default: throw "SymmetricMatrix::Inverse error: Unsupported matrix size";
        }
    }
    
    
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
    
}

#endif
