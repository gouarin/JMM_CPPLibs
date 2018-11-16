//
//  SymmetricMatrixType.hxx
//  
//
//  Created by Jean-Marie Mirebeau on 09/02/2018.
//

#ifndef SymmetricMatrixType_hxx
#define SymmetricMatrixType_hxx

template<typename TC, size_t VD>  int SymmetricMatrix<TC,VD>::
LinearizedIndex(int i, int j){
    assert(IsInRange(i, j));
    if(i<j) std::swap(i, j);
    return (i*(i+1))/2+j;
}

template<typename TC, size_t VD> template<typename T> auto SymmetricMatrix<TC,VD>::
ScalarProduct(const Vector<T, Dimension> & u, const Vector<T, Dimension> & v) const
-> ComponentType
{
    ComponentType sum=0;
    for(int i=0; i<Dimension; ++i){
        ComponentType sumi=0;
        for(int j=0; j<Dimension; ++j)
            sumi+=coef(i,j)*ComponentType(v[j]);
            sum+=ComponentType(u[i])*sumi;
            }
    return sum;
}

template<typename TC, size_t VD> template<typename T> auto SymmetricMatrix<TC,VD>::
operator*(const Vector<T,Dimension> & u) const -> VectorType
{
    VectorType v;
    for(int i=0; i<Dimension; ++i){
        v[i]=0;
        for(int j=0; j<Dimension; ++j)
            v[i]+=coef(i, j)*ComponentType(u[j]);
            }
    return v;

}

template<typename TC, size_t VD> auto SymmetricMatrix<TC,VD>::
Zero() -> SymmetricMatrix
{
    SymmetricMatrix m;
    std::fill(m.data.begin(), m.data.end(), ComponentType(0));
    return m;
}

template<typename TC, size_t VD> auto SymmetricMatrix<TC,VD>::
Identity() -> SymmetricMatrix
{
    SymmetricMatrix m;
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<=i; ++j)
            m(i,j)=(i==j);
    return m;
}

template<typename TC, size_t VD> auto SymmetricMatrix<TC,VD>::
Diagonal(const VectorType & u) -> SymmetricMatrix
{
    SymmetricMatrix m;
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<=i; ++j)
            m(i,j)= i==j ? u[i] : ComponentType(0);
    return m;
}

template<typename TC, size_t VD> template<typename T> auto SymmetricMatrix<TC,VD>::
RankOneTensor(const Vector<T,Dimension> & u) -> SymmetricMatrix {
    SymmetricMatrix m;
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<=i; ++j)
            m(i,j)=u[i]*u[j];
    return m;
}

template<typename TC,size_t VD> auto SymmetricMatrix<TC,VD>::
RandomPositive() -> SymmetricMatrix {    
    MatrixType a=MatrixType::Random();
    return FromUpperTriangle(a.Transpose()*a);
}

template<typename TC,size_t VD> auto SymmetricMatrix<TC,VD>::
FromUpperTriangle(const MatrixType & mat) -> SymmetricMatrix
{
    SymmetricMatrix m;
    for(int i=0; i<Dimension; ++i)
        for(int j=i; j<Dimension; ++j)
            m(i,j)=mat(i,j);
            return m;
}

template<typename TC,size_t VD> auto SymmetricMatrix<TC,VD>::
operator=(ComponentType a) -> SymmetricMatrix
{
    static_assert(Dimension==1,"Assignment from scalar to matrix in dimension 1 only");
    data[0]=a;
    return *this;
}


template<typename TC, size_t VD> template<typename T, size_t D> auto SymmetricMatrix<TC,VD>::
Gram(const std::array<Vector<T,Dimension>, D> & a) const
-> SymmetricMatrix<ComponentType,D> {
    SymmetricMatrix<ComponentType,D> m;
    for(int i=0; i<D; ++i){
        const Vector<ComponentType, Dimension> mai = operator*(a[i]);
        for(int j=0; j<=i; ++j)
            m(i,j) = mai.ScalarProduct(a[j]);
    }
    /*
     SymmetricMatrix<ComponentType,D> m;
     for(int i=0; i<D; ++i)
     for(int j=0; j<=i; ++j)
     m(i,j) = ScalarProduct(a[i],a[j]);*/
    return m;
}

template<typename TC, size_t VD> template<size_t D> auto SymmetricMatrix<TC,VD>::
Gram(const Matrix<ComponentType,Dimension,D> & a) const
-> SymmetricMatrix<ComponentType, D>
{
    Matrix<ComponentType,D,Dimension> ta = a.Transpose();
    return Gram(*(const std::array<Vector<ComponentType,Dimension>,D>*)&ta);
}

template<typename TC, size_t VD> template<size_t D> auto SymmetricMatrix<TC,VD>::
GramT(const Matrix<ComponentType,D,Dimension> & a) const
-> SymmetricMatrix<ComponentType, D>
{
    return Gram(*(const std::array<Vector<ComponentType,Dimension>,D>*)&a);
}

template<typename TC, size_t VD> template<typename T, size_t D> auto SymmetricMatrix<TC,VD>::
EuclideanGram(const std::array<Vector<T,D>, Dimension> & a) -> SymmetricMatrix
{
    SymmetricMatrix m;
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<=i; ++j)
            m(i,j) = a[i].ScalarProduct(a[j]);
    return m;
}

// Conversion
template<typename TC,size_t VD> SymmetricMatrix<TC,VD>::
operator MatrixType() const {
    MatrixType result;
    for(int i=0; i<Dimension; ++i)
        for(int j=0; j<Dimension; ++j)
            result(i,j)=coef(i, j);
    return result;
}

template<typename TC,size_t VD> template<typename T> auto SymmetricMatrix<TC,VD>::
CastCoordinates(const SymmetricMatrix<T, Dimension> & m0) -> SymmetricMatrix {
    SymmetricMatrix m;
    m.data = Vector<ComponentType,InternalDimension>::CastCoordinates(m0.data);
    return m;
}


// Matrix operations
template<typename TC,size_t VD> auto SymmetricMatrix<TC,VD>::
Trace() const -> ComponentType
{
    ComponentType result=0;
    for(int i=0; i<Dimension; ++i)
        result+=coef(i,i);
        return result;
}


template<typename TC, size_t VD>
TC
SymmetricMatrix<TC, VD>::Determinant() const {
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
        default: {
            LinearAlgebra::Matrix<TC, VD, VD> mat;
            for(int i=0; i<VD; ++i)
                for(int j=0; j<VD; ++j)
                    mat(i,j)=coef(i,j);
            return mat.Determinant();
        }
    }
}

template<typename TC, size_t VD>
SymmetricMatrix<TC, VD>
SymmetricMatrix<TC, VD>::Inverse() const {
    if(VD>3) return FromUpperTriangle(this->operator MatrixType().Inverse());
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
                    const int i1 = (i+1)%3, i2=(i+2)%3, j1=(j+1)%3, j2=(j+2)%3;
                    m(i,j) = coef(i1,j1)*coef(i2,j2)-coef(i1,j2)*coef(i2,j1);
                }
            m/=det; return m;
        default:
            assert(false);
			ExceptionMacro("SymmetricMatrix::Inverse() error : unsupported dimension");
    }
}

template<typename TC, size_t VD> auto SymmetricMatrix<TC,VD>::
CGSolve(VectorType r) const -> VectorType {
    VectorType x= VectorType::Constant(0);
    const SymmetricMatrix & m = *this;
    VectorType p = r, r0, q;
    ComponentType r2=r.SquaredNorm(), r02=0; //Dummy for r02.
    for(int k=0; k<Dimension && r2>0; ++k){
        if(k>0){
            const ComponentType beta = r2/r02;
            p=r+beta*p;
        }
        q = m*p;
        const ComponentType pq = p.ScalarProduct(q);
        if(pq == 0) return x;
        const ComponentType alpha = r.SquaredNorm()/pq;
        x+=alpha*p;
        r0=r;
        r02=r2;
        r-=alpha*q;
        r2=r.SquaredNorm();
    }
    return x;
}

#endif /* SymmetricMatrixType_hxx */
