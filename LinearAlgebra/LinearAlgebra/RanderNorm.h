//
//  RanderNorm.h
//  LinearAlgebra
//
//  Created by Jean-Marie Mirebeau on 16/04/2016.
//
//

#ifndef RanderNorm_LinearAlgebra_h
#define RanderNorm_LinearAlgebra_h

#include "VectorType.h"
#include "SymmetricMatrixType.h"

namespace LinearAlgebra {
    
template<typename TScalar, size_t VDimension>
struct RanderNorm {
    
    typedef TScalar ScalarType;
    static const size_t Dimension = VDimension;

    typedef Vector<ScalarType,Dimension> VectorType;
    typedef SymmetricMatrix<ScalarType,Dimension> SymmetricMatrixType;
    
    
    SymmetricMatrixType m;
    VectorType w;
    
    ScalarType Norm(const VectorType & v) const {return m.Norm(v) - w.ScalarProduct(v);}
    const VectorType Gradient(const VectorType & v) const {return m*v/m.Norm(v) - w;}
    bool IsAcute(const VectorType & u, const VectorType & v) const {
        return Gradient(u).ScalarProduct(v)>=0 && Gradient(v).ScalarProduct(u)>=0;}
    
    const RanderNorm DualNorm() const;
    bool IsDefinite() const {return m.Inverse().SquaredNorm(w)<1;}
};
    
template<typename TS, size_t VD>
const RanderNorm<TS,VD> RanderNorm<TS,VD>::DualNorm() const {
    const SymmetricMatrixType s = (m-SymmetricMatrixType::RankOneTensor(w)).Inverse();
    const VectorType omega = s*w;
    return RanderNorm{s*(1+w.ScalarProduct(omega)),omega};
}

template<typename TS, size_t VD>
std::ostream & operator << (std::ostream & os, const RanderNorm<TS,VD> & norm){
    return os << "{" << norm.m << "," << norm.w << "}";
}
    
}

#endif
