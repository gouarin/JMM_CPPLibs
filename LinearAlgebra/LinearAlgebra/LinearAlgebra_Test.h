//
//  LinearAlgebra_Test.h
//  AmongConvex2
//
//  Created by Jean-Marie Mirebeau on 24/07/13.
//  Copyright (c) 2013 Jean-Marie Mirebeau. All rights reserved.
//

#ifndef AmongConvex2_LinearAlgebra_Test_h
#define AmongConvex2_LinearAlgebra_Test_h

#include "AffineTransformType.h"
bool LinearAlgebra_Test(){
    using std::cout;
    using std::endl;

    {
        typedef double ComponentType;
        typedef LinearAlgebra::AffineTransform<int, 3, 2> AffineTransformType;
        
        AffineTransformType::InputPointType P2{0,1};
        AffineTransformType::InputVectorType V2{1,2};
        
        AffineTransformType::OutputPointType P3{-1,0,1};
        AffineTransformType::OutputVectorType V3{-1,0,-2};
        
        
        cout << P2 << P2.AreAllCoordinatesNonNegative() << P2.AreAllCoordinatesPositive() << endl;
        cout << V2 << V2.AreAllCoordinatesNonNegative() << V2.AreAllCoordinatesPositive() << endl;
        cout << P3 << P3.AreAllCoordinatesNonNegative() << P3.AreAllCoordinatesPositive() << endl;
        cout << V3 << V3.AreAllCoordinatesNonNegative() << V3.AreAllCoordinatesPositive() << endl;
        
        cout << P2+V2 << V2+V2 << P2-P2 << endl;
        cout << P3+V3 << V3+V3 << P3-V3 << endl;
        
        AffineTransformType A;
        A.linearPart = AffineTransformType::MatrixType::FromRows({{ {{1,2}}, {{3,4}}, {{5,6}} }});
        A.imageOfOrigin = P3;
        
        cout << A(P2) << A.linearPart*V2 << endl;
        cout << A.linearPart << endl;
    }
    
    {

        typedef LinearAlgebra::AffineTransform<double,2,2> AffineTransformType;
        typedef AffineTransformType::MatrixType MatrixType;
        typedef AffineTransformType::OutputPointType PointType;
        
        AffineTransformType A(MatrixType::Rotation(1), PointType{1,2});

        cout << "A: " << A << endl;
        
        AffineTransformType B=A.Inverse();
        
        cout << "B: " << B << endl;
        cout << "R(-1): "<< MatrixType::Rotation(-1)<<endl;
        
        cout << A.Compose(B) << " " << B.Compose(A) << endl;
    }
    return true;
}

#endif
