//
//  itkImage.h
//  MatlabInterface
//
//  Created by Jean-Marie Mirebeau on 31/03/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//


// !! THIS IS NOT ITK !!
// This is a quick mockup made to reproduce only very basic functionality, needed for my C++ to Matlab interface.

#ifndef MatlabInterface_itkImage_h
#define MatlabInterface_itkImage_h

#include <array>
#include <vector>
#include <numeric>
#include <map>
#include "MatlabITKInterface/MexMessageWrapper.h"

// **************** Exceptions ******************


#define itkExceptionMacro(x) \
{ \
std::ostringstream message; \
message << "File: " << __FILE__ << ", Line: " << __LINE__ << ", PseudoItk::ERROR: " << this->GetNameOfClass() << "(" << this << "): " x; \
throw std::logic_error(message.str()); \
}

#define itkGenericExceptionMacro(x)  \
{ \
std::ostringstream message;  \
message << "File: " << __FILE__ << ", Line: " << __LINE__ << "PseudoItk::ERROR: " x; \
throw std::logic_error(message.str()); \
}

namespace itk {
    
    // ******************** Reference counted object ********************
    template<typename Object>
    struct Pointer {
        operator Object*(){return m_Pointer;}
        Object * operator->() const {return m_Pointer;}
        //Object & operator*() const {return *m_Pointer;}
        void operator=(Object*pointer){
            if(pointer==m_Pointer) return;
            if(m_Pointer!=nullptr){
                auto it = referenceCount.find(m_Pointer);
                if(it==referenceCount.end()) {itkExceptionMacro(<<"Pointer referenceCount error : no pointer.");}
                else if(it->second<=0) {itkExceptionMacro(<<"Pointer referenceCount error : count is zero.");}
                else if(it->second==1) {referenceCount.erase(it); delete m_Pointer;} //std::cout << "Deleting" << m_Pointer <<"\n";
                else it->second = it->second-1;
            }
            if(pointer!=nullptr){
                auto it = referenceCount.find(pointer);
                if(it==referenceCount.end()) referenceCount.insert({pointer,1});
                else it->second = it->second+1;
            }
            m_Pointer=pointer;
        }
        void operator=(const Pointer & ptr){operator=(ptr.m_Pointer);}
        Pointer(){};
        Pointer(const Pointer & ptr){operator=(ptr);};
        ~Pointer(){operator=(nullptr);}
        std::string GetNameOfClass() const {return "Pointer (Reference counted)";}
        bool IsNull() const {return m_Pointer==nullptr;}
        void PrintSelf(std::ostream & os) const {os << "Pointer("; if(m_Pointer!=nullptr) os<<*m_Pointer; os << ")";}
        friend std::ostream & operator <<(std::ostream & os, Pointer x){ x.PrintSelf(os); return os;}
        //std::ostream & operator << (std::ostream & os) const {
        //    os << "Pointer(" << (m_Pointer!=nullptr ? *m_Pointer : "nullptr") << ")"; return os;}
    protected:
        Object * m_Pointer = nullptr;
        static std::map<Object*,int> referenceCount;
    };
    template<typename Object> std::map<Object*,int> Pointer<Object>::referenceCount;
//    template<typename Object> std::ostream & operator << (std::ostream & os, Pointer<Object>) {return os;};
    // ****************** Basic types ****************
    
    typedef long DiscreteType;
    typedef double ScalarType;

    template<typename TComponent,int NDimension>
    struct Point : std::array<TComponent, NDimension> {
        typedef TComponent ComponentType;
        static const int Dimension = NDimension;
        void Fill(ComponentType c){this->fill(c);}
        Point(ComponentType c){Fill(c);}
        Point(){};
        friend std::ostream & operator << (std::ostream & os, Point x) {
            os << "{"; for(int i=0; i<Dimension; ++i) os << x[i] << ","; os << "}"; return os;}
    };
    
    // ****************** Images ******************
    
    template<DiscreteType NDimension>
    struct ImageRegion {
        static const DiscreteType ImageDimension = NDimension;
        typedef Point<DiscreteType,ImageDimension> IndexType;
        typedef IndexType SizeType;
        void SetIndex(IndexType index){m_Index=index;}
        void SetSize(SizeType size){m_Size=size;}
        const IndexType & GetIndex() const {return m_Index;}
        const SizeType & GetSize() const {return m_Size;}
        std::string GetNameOfClass() const {return "ImageRegion";}
        ImageRegion(){};
        ImageRegion(IndexType index, SizeType size):m_Index(index),m_Size(size){};
        bool operator==(const ImageRegion & other) const {return m_Index==other.m_Index && m_Size==other.m_Size;}
        bool operator!=(const ImageRegion & other) const {return !operator==(other);}
        friend std::ostream & operator << (std::ostream & os, ImageRegion x) {os << "Region(" << x.m_Size << ")"; return os;}
    protected:
        IndexType m_Index; SizeType m_Size;
    };
    
    template<typename TPixel, DiscreteType NDimension>
    struct Image {
        typedef TPixel PixelType;
        static const DiscreteType ImageDimension = NDimension;
        typedef ImageRegion<ImageDimension> RegionType;
        typedef typename RegionType::IndexType IndexType;
        typedef typename RegionType::SizeType SizeType;
        typedef Point<ScalarType,ImageDimension> PointType;
        typedef PointType SpacingType;
        typedef Pointer<Image> Pointer;
        
        void SetRegions(RegionType region){m_Region=region;}
        const RegionType & GetBufferedRegion() const {return m_Region;}
        const RegionType & GetRequestedRegion() const {return m_Region;}
        const RegionType & GetLargestPossibleRegion() const {return m_Region;}
        
        void SetSpacing(SpacingType spacing){m_Spacing=spacing;}
        const SpacingType & GetSpacing() const {return m_Spacing;}
        
        PixelType & GetPixel(const IndexType & index) {return buffer[LinearIndex(index)];}
        void SetOrigin(PointType origin){m_Origin=origin;}
        const PointType GetOrigin() const {return m_Origin;}
        
        static Pointer New() {Image * raw = new Image; Pointer pointer; pointer=raw; return pointer;}
        void Allocate(){buffer.resize(NumberOfPixels());}
              PixelType * GetBufferPointer()       {return &buffer[0];}
        const PixelType * GetBufferPointer() const {return &buffer[0];}
              std::vector<PixelType> * GetPixelContainer()       {return &buffer;}
        const std::vector<PixelType> * GetPixelContainer() const {return &buffer;}
        std::string GetNameOfClass() const {return "Image";}
        friend std::ostream & operator << (std::ostream & os, Image x) {os << "Image(" << x.m_Region << ",bufferSize"<<x.buffer.size()<<")"; return os;}
    protected:
        RegionType m_Region;
        SpacingType m_Spacing;
        PointType m_Origin;
        std::vector<PixelType> buffer;
        DiscreteType NumberOfPixels(){
            return std::accumulate(m_Region.GetSize().begin(), m_Region.GetSize().end(),DiscreteType(1),
                                   [](DiscreteType a, DiscreteType b){return a*b;});}
        DiscreteType LinearIndex(const IndexType & index) const {
            const auto size = m_Region.GetSize();
            const auto rIndex = m_Region.GetIndex();
            DiscreteType result=0;
            for(int i=ImageDimension-1; i>=0; --i){
                const DiscreteType indexI = index[i]-rIndex[i];
                assert(indexI>=0 && indexI<size[i]);
                if(i!=ImageDimension-1) result*=size[i];
                result+=indexI;}
            return result;}
        Image(){};
        template<typename TImage> friend class ImageRegionConstIteratorWithIndex;
    };
}


#endif
