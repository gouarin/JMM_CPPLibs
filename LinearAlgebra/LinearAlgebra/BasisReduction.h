// Copyright 2017 Jean-Marie Mirebeau, University Paris-Sud, CNRS, University Paris-Saclay
// Distributed WITHOUT ANY WARRANTY. Licensed under the Apache License, Version 2.0, see http://www.apache.org/licenses/LICENSE-2.0

#ifndef BasisReduction_h
#define BasisReduction_h

#include "SymmetricMatrixType.h"

namespace LinearAlgebra {
    
template<typename TScalar, typename TDiscrete, size_t VDimension>
struct BasisReduction {
    typedef TScalar ScalarType;
    typedef TDiscrete DiscreteType;
    static const size_t Dimension = VDimension;
    typedef Vector<ScalarType, Dimension> VectorType;
    typedef Vector<DiscreteType, Dimension> DiscreteVectorType;
    typedef SymmetricMatrix<ScalarType, Dimension> SymmetricMatrixType;
    
    typedef std::array<DiscreteVectorType, Dimension>   BasisType;
    typedef std::array<DiscreteVectorType, Dimension+1> SuperbaseType;
    typedef std::array<std::pair<DiscreteVectorType, ScalarType>,
    (Dimension*(Dimension+1))/2> TensorDecompositionType;
    
    static const BasisType CanonicalBasis();
    static const SuperbaseType CanonicalSuperBase();
    static void ReducedBasis(const SymmetricMatrixType &, BasisType &);
    static void ObtuseSuperbase(const SymmetricMatrixType &, SuperbaseType &);
    static TensorDecompositionType TensorDecomposition(const SymmetricMatrixType &);
    static int maxIt;
protected:
    template<size_t VD, typename Dummy=void> struct TensorDecompositionHelper;
};
template<typename TS,typename TD,size_t VD> int BasisReduction<TS,TD,VD>::maxIt=200;
#include "BasisReduction.hpp"
}


#endif /* BasisReduction_h */
