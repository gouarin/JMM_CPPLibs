// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef AsymmetricQuadraticNorm_h
#define AsymmetricQuadraticNorm_h

/*
 
 Implements the norm v.M.v+max(w.v,0)^2
 
 */

#include "SymmetricMatrixType.h"

namespace LinearAlgebra {
    
template<typename TScalar, size_t VDimension>
struct AsymmetricQuadraticNorm {
    
    typedef TScalar ScalarType;
    static const size_t Dimension = VDimension;
    
    typedef Vector<ScalarType, Dimension> VectorType;
    typedef SymmetricMatrix<ScalarType, Dimension> SymmetriMatrixType;
    
    SymmetriMatrixType m;
    VectorType w;
    
    ScalarType Norm(const VectorType & v) const {
        const ScalarType scal = ScalPos(v);
        return sqrt(m.SquaredNorm(v)+scal*scal);
    }
    
     //Positive multiple of gradient
    VectorType MGrad(const VectorType & v) const {return m*v+ScalPos(v)*w;}
    VectorType Grad(const VectorType & v) const {
        const VectorType g = MGrad(v);
        const ScalarType squaredNorm = g.ScalarProduct(v);
        return g/sqrt(squaredNorm);
    }
    
    // Equivalent to u.Grad(v)>=0 && v.Grad(u)>=0.
    bool IsAcute(const VectorType & u, const VectorType & v) const {
        const ScalarType
        muv = m.ScalarProduct(u,v), wu = w.ScalarProduct(v), wv=w.ScalarProduct(v);
        return muv+std::min(wu*std::max(wv,0.),wv*std::max(wu,0.)) >= 0.;
    }
#pragma todo("Implement dual")
    const AsymmetricQuadraticNorm DualNorm() const {throw "AQN dual is TODO";}
    bool IsDefinite() const {return m.IsDefinite();}
    
protected:
    ScalarType ScalPos(const ScalarType & v) const {return std::max(w.ScalarProduct(v),0);}
};
}

#endif /* AsymmetricQuadraticNorm_h */
