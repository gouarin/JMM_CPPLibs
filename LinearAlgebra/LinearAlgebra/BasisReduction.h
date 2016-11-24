//
//  BasisReduction.h
//  ExternLFM
//
//  Created by Jean-Marie Mirebeau on 21/10/2016.
//
//

#ifndef BasisReduction_h
#define BasisReduction_h

#include "SymmetricMatrixType.h"

namespace LinearAlgebra {
    
template<typename TScalar, typename TDiscrete, size_t VDimension>
struct BasisReduction {
    static_assert(VDimension<=3,"Sorry, dimensions >=4 are not supported");
    typedef TScalar ScalarType;
    typedef TDiscrete DiscreteType;
    static const size_t Dimension = VDimension;
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
protected:
    template<size_t VD, typename Dummy=void> struct TensorDecompositionHelper;
};
#include "BasisReduction.hpp"
}


#endif /* BasisReduction_h */
