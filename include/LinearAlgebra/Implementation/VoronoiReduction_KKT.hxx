//
//  VoronoiReduction_KKT.hxx
//  Voronoi45
//
//  Created by Jean-Marie Mirebeau on 06/02/2018.
//  Copyright © 2018 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef VoronoiReduction_KKT_h
#define VoronoiReduction_KKT_h

template<typename TS, int VD> auto VoronoiFirstReduction<TS,VD>::
TensorDecomposition(const SymmetricMatrixType & m,ScalarType tol) -> KKTRelationType {
    SimplexStateType state(m);
    if(tol>=0){GreedyBasis(state,tol);}
    Minimize(state);
    return KKT(state);
}


template<typename TS, int VD> auto VoronoiFirstReduction<TS,VD>::
KKT(const SimplexStateType & state) -> KKTRelationType {
    KKTRelationType result;
    typedef LinearAlgebra::Matrix<ScalarType, SymDimension, SymDimension> KKTMat;
    typedef typename KKTMat::InputVectorType KKTVec;
    if constexpr(Dimension==1){
        result.offsets[0][0]=1;
        result.weights[0]=state.m(0,0);
    }
    if constexpr(Dimension==2){
        constexpr KKTMat A{1,1,0,0,1,1,0,-1,0};
        constexpr std::array<VectorType,SymDimension> support = {{
            {1,0},{0,1},{1,-1}
        }};
        result.weights = A*KKTVec(state.m.data);
        const MatrixType aInv = state.a.Inverse();
        for(int i=0; i<SymDimension; ++i){
            result.offsets[i] = OffsetType::CastCoordinates(aInv*support[i]);}
    }
    if constexpr(Dimension==3){
        constexpr KKTMat A{
            1,1,0,1,0,0,0,1,1,0,1,0,0,0,0,1,1,1,0,0,0,-1,
                0,0,0,-1,0,0,0,0,0,0,0,0,-1,0
        };
        constexpr std::array<VectorType,SymDimension> support = {{
            {1,0,0},{0,1,0},{0,0,1},{1,0,-1},{1,-1,0},{0,1,-1}
        }};
        result.weights = A*KKTVec(state.m.data);
        const MatrixType aInv = state.a.Inverse();
        for(int i=0; i<SymDimension; ++i){
            result.offsets[i] = OffsetType::CastCoordinates(aInv*support[i]);}
    }
    
    if constexpr(Dimension==4){
        if(state.vertex==1){
            constexpr KKTMat A = {1,1,0,1,0,0,1,0,0,0,0,1,1,0,1,0,0,1,0,0,0,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,-1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0};
            constexpr std::array<VectorType, SymDimension> support = {{
                {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,0,0,-1},
                {1,0,-1,0},{1,-1,0,0},{0,1,0,-1},{0,1,-1,0},{0,0,1,-1}
            }};
            result.weights = A*KKTVec(state.m.data);
            const MatrixType aInv = state.a.Inverse();
            for(int i=0; i<SymDimension; ++i){
                result.offsets[i] = OffsetType::CastCoordinates(aInv*support[i]);}
        } else { assert(state.vertex==0);
            constexpr KKTMat A = {1,1,0,1,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,0
                ,0,1,1,1,1,0,1,0,0,0,0,0,0,0,0,1,1,1,0,0,0,-1,0,0,-1,0,0,0,0,
                -1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,-1,0,0,0,0,0,0
                ,0,0,0,0,0,-1,0,-1,0,0,0,0,0,0,0,1,0,0,0};
            constexpr std::array<VectorType, 12> support = {{
            {1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1},{1,0,-1,0},{1,-1,0,0},
                {0,1,0,-1},{0,1,-1,0},{0,0,1,-1},{1,0,-1,1},{1,-1,0,1},{1,-1,-1,1}
            }};
            
            // Compute a particular solution to the linear system, with possibly negative weights
            const KKTVec sol0 = A*KKTVec(state.m.data);
            const ScalarType maxSol0 = *std::max_element(sol0.begin(), sol0.end(),
                [](ScalarType a, ScalarType b)->ScalarType{
                return std::max<ScalarType>(a,std::abs(b));});
            
            // Call to a linear solver, to find positive weights
            const int d = 2;
            const int max_size = 13;
            const int m = max_size;
            
            ScalarType halves[max_size][d+1] = { // TODO : solve by hand
                // Ask that the ten first multipliers be positive
                {0,1,0},
                {1,0,0},
                {-1,-1,0},
                {0,1,0},
                {1,0,0},
                {-1,-1,0},
                {-1,-1,0},
                {0,1,0},
                {1,0,0},
                {-1,-1,0},
                
                // last two must be positive as well
                {1,0,0},
                {0,1,0},
                
                // projective component positive
                {0,0,1}
            };
            for(int i=0; i<SymDimension; ++i)
                halves[i][2] = sol0[i]/maxSol0;
            
            // Minimize the sum of all coefficients (check)
            ScalarType n_vec[d+1] = {1,1,0};
            ScalarType d_vec[d+1] = {0,0,1};
            
            ScalarType opt[d+1];
            ScalarType work[((max_size+3)*(d+2)*(d-1))/2];
            
            const int BadIndex = 1234567890;
            int next[max_size] = {1,2,3,4,5,6,7,8,9,10,11,12,13};
            int prev[max_size] = {BadIndex,0,1,2,3,4,5,6,7,8,9,10,11};
            
            dlinprog(&halves[0][0], 0, m, n_vec, d_vec, d, opt, work, next, prev, max_size);
            
            // TODO : check that status is correct
            // Get the solution, and find the non-zero weights, which should be positive.
            const std::array<ScalarType, 2> coef{{opt[0]/opt[2],opt[1]/opt[2]}};
            std::array<ScalarType, 12> weights;
            for(int i=0; i<SymDimension; ++i){
                weights[i] = maxSol0*(
                coef[0]*halves[i][0]+coef[1]*halves[i][1]+halves[i][2]);
            }
            weights[10] = maxSol0*coef[0];
            weights[11] = maxSol0*coef[1];
            
            std::array<int, 12> ord = {{0,1,2,3,4,5,6,7,8,9,10,11}};
            std::nth_element(ord.begin(), ord.begin()+2, ord.end(),
                [&weights](int i, int j)->bool {return weights[i]<weights[j];});
            
            const MatrixType aInv = state.a.Inverse();
            for(int i=0; i<SymDimension; ++i){
                const int j=ord[i+2];
                result.weights[i] = weights[j];
                result.offsets[i] = OffsetType::CastCoordinates(aInv*support[j]);
            } // for i
        } // if state.vertex == 0
    } // dimension == 4
                
    if constexpr(Dimension==5){
        if(state.vertex==1){
            constexpr KKTMat A = {1,1,0,1,0,0,1,0,0,0,1,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0};
            constexpr std::array<VectorType, SymDimension> support = {{
                {1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1},{1,0,0,0,-1},{1,0,0,-1,0},{1,0,-1,0,0},{1,-1,0,0,0},{0,1,0,0,-1},{0,1,0,-1,0},{0,1,-1,0,0},{0,0,1,0,-1},{0,0,1,-1,0},{0,0,0,1,-1}
            }};
            result.weights = A*KKTVec(state.m.data);
            const MatrixType aInv = state.a.Inverse();
            for(int i=0; i<SymDimension; ++i){
                result.offsets[i] = OffsetType::CastCoordinates(aInv*support[i]);}
        } else if(state.vertex==2){
            constexpr KKTMat A = {1,0.5,0,0.5,-0.5,0,-1,0,0,0,-1,0,0,0.5,0,0,0.5,1,-0.5,0.5,0,0,-1,0,0,0,-1,0,0.5,0,0,-0.5,0,0.5,0.5,1,0,0,-1,0,0,0,-1,0.5,0,0,0.5,0,0.5,0.5,0,-1,-1,-1,1,0,0,0,0.5,0,0,0.5,0,0.5,0.5,0,0,0,0,0,-1,-1,-1,0.5,1,0,-0.5,0,-0.5,0.5,0,1,0,0,0,0,0,0,-0.5,0,0,-0.5,0,-0.5,0.5,0,0,0,0,0,1,0,0,-0.5,0,0,-0.5,0,0.5,-0.5,0,0,1,0,0,0,0,0,-0.5,0,0,-0.5,0,0.5,-0.5,0,0,0,0,0,0,1,0,-0.5,0,0,0.5,0,-0.5,-0.5,0,0,0,1,0,0,0,0,-0.5,0,0,0.5,0,-0.5,-0.5,0,0,0,0,0,0,0,1,-0.5,0,0,0.5,0,-0.5,-0.5,0,0,0,0,0,0,0,0,0.5,0,0,-0.5,0,0.5,-0.5,0,0,0,0,0,0,0,0,0.5,0,0,-0.5,0,-0.5,0.5,0,0,0,0,0,0,0,0,0.5,0,0,0.5,0,0.5,0.5,0,0,0,0,0,0,0,0,-0.5,0};
            constexpr std::array<VectorType, SymDimension> support = {{
                {1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1},{1,0,0,1,0},{1,0,0,0,1},{0,1,0,1,0},{0,1,0,0,1},{0,0,1,1,0},{0,0,1,0,1},{1,1,0,1,1},{1,0,1,1,1},{0,1,1,1,1},{1,1,1,1,1}
                }};
            result.weights = A*KKTVec(state.m.data);
            const MatrixType aInv = state.a.Inverse();
            for(int i=0; i<SymDimension; ++i){
                result.offsets[i] = OffsetType::CastCoordinates(aInv*support[i]);}
        } else { assert(state.vertex==0);
            constexpr KKTMat A = {1,1,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,-1,0,0,0,-1,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0};
            constexpr std::array<VectorType, 20> support = {{
                {1,0,0,0,0},{0,1,0,0,0},{0,0,1,0,0},{0,0,0,1,0},{0,0,0,0,1},{1,0,0,-1,0},{1,0,-1,0,0},{1,-1,0,0,0},{0,1,0,0,-1},{0,1,0,-1,0},{0,1,-1,0,0},{0,0,1,0,-1},{0,0,1,-1,0},{0,0,0,1,-1},{1,0,0,-1,1},{1,0,-1,0,1},{1,-1,0,0,1},{1,0,-1,-1,1},{1,-1,0,-1,1},{1,-1,-1,0,1}
            }};
            
            // Compute a particular solution to the linear system, with possibly negative weights
            const KKTVec sol0 = A*KKTVec(state.m.data);
            const ScalarType maxSol0 = *std::max_element(sol0.begin(), sol0.end(),
                                                         [](ScalarType a, ScalarType b)->ScalarType{
                                                             return std::max<ScalarType>(a,std::abs(b));});
            
            // Call to a linear solver, to find positive weights
            const int d = 5;
            const int max_size = 21;
            const int m = max_size;
            
            // TODO : many constraints appear several times.
            // (Only 5 distinct and non axis related.) Simplify ?
            ScalarType halves[max_size][d+1] = {
                // Ask that the fifteen first multipliers be positive
                {0,0,1,1,1,0},
                {0,1,0,0,0,0},
                {1,0,0,0,0,0},
                {-1,-1,-1,-1,-1,0},
                {0,0,1,1,1,0},
                {1,1,0,0,1,0},
                {-1,0,-1,0,-1,0},
                {0,-1,0,-1,-1,0},
                {0,-1,0,-1,-1,0},
                {0,0,0,1,0,0},
                {0,0,0,0,1,0},
                {-1,0,-1,0,-1,0},
                {0,0,1,0,0,0},
                {1,1,0,0,1,0},
                {-1,-1,-1,-1,-1,0},
                
                // last five must be positive as well
                {1,0,0,0,0,0},
                {0,1,0,0,0,0},
                {0,0,1,0,0,0},
                {0,0,0,1,0,0},
                {0,0,0,0,1,0},
                // projective component positive
                {0,0,0,0,0,1}
            };
            for(int i=0; i<SymDimension; ++i)
                halves[i][5] = sol0[i]/maxSol0;
                
                // Minimize the sum of all coefficients (check)
                ScalarType n_vec[d+1] = {1,1,1,1,1,0};
                ScalarType d_vec[d+1] = {0,0,0,0,0,1};
                
                ScalarType opt[d+1];
            ScalarType work[((max_size+3)*(d+2)*(d-1))/2];
            
            const int BadIndex = 1234567890;
            int next[max_size] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21};
            int prev[max_size] = {BadIndex,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
            
            dlinprog(&halves[0][0], 0, m, n_vec, d_vec, d, opt, work, next, prev, max_size);
            
            // TODO : check that status is correct
            // Get the solution, and find the non-zero weights, which should be positive.
            const std::array<ScalarType, 5> coef{{
                opt[0]/opt[5],opt[1]/opt[5],opt[2]/opt[5],opt[3]/opt[5],opt[4]/opt[5]}};
            std::array<ScalarType, 20> weights;
            for(int i=0; i<SymDimension; ++i){
                ScalarType & w = weights[i];
                w=0;
                for(int j=0; j<5; ++j) {w+=coef[j]*halves[i][j];}
                w+=halves[i][5];
                w*=maxSol0;
            }
            for(int i=0; i<5; ++i)
                weights[15+i] = maxSol0*coef[i];
            
                std::array<int, 20> ord = {{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}};
            std::nth_element(ord.begin(), ord.begin()+5, ord.end(),
                             [&weights](int i, int j)->bool {return weights[i]<weights[j];});
            
            const MatrixType aInv = state.a.Inverse();
            for(int i=0; i<SymDimension; ++i){
                const int j=ord[i+5];
                result.weights[i] = weights[j];
                result.offsets[i] = OffsetType::CastCoordinates(aInv*support[j]);
            } // for i
        }
    }
    return result;
}

template<typename TS, int VD> void
VoronoiFirstReduction<TS, VD>::KKTRelationType::
PrintSelf(std::ostream & os) const {
    os << "{"
    ExportArrayArrow(weights)
    ExportArrayArrow(offsets)
    << "}";
}

#endif /* VoronoiReduction_KKT_h */
