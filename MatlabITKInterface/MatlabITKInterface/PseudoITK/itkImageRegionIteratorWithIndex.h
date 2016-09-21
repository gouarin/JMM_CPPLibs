//
//  itkImageRegionIteratorWithIndex.h
//  MatlabInterface
//
//  Created by Jean-Marie Mirebeau on 31/03/2015.
//  Copyright (c) 2015 Jean-Marie Mirebeau. All rights reserved.
//

// !! THIS IS NOT ITK !!
// This is a quick mockup made to reproduce only very basic functionality, needed for my C++ to Matlab interface.

#ifndef MatlabInterface_itkImageRegionIteratorWithIndex_h
#define MatlabInterface_itkImageRegionIteratorWithIndex_h

#include "itkImage.h"

namespace itk {

template<typename TImage>
struct ImageRegionConstIterator {
    typedef TImage ImageType;
    typedef typename ImageType::Pointer ImagePointer;
    typedef typename ImageType::RegionType ImageRegionType;
    typedef typename ImageType::PixelType PixelType;

    ImageRegionConstIterator(const ImageType * ptr,const ImageRegionType & region):m_Ptr(ptr){
        if(ptr->GetBufferedRegion()!=region)
            itkExceptionMacro( << "ImageRegionConstIterator error: region mismatch");
        end=ptr->GetBufferPointer()+ptr->GetPixelContainer()->size();
    }
    const PixelType & Value() const {return *current;}
    void GoToBegin(){current=m_Ptr->GetBufferPointer();}
    bool IsAtEnd(){return current==end;}
    void operator++(){++current;}
    std::string GetNameOfClass() const {return "ImageRegionConstIterator";}
protected:
    const ImageType * m_Ptr=nullptr;
    const PixelType * end=nullptr, * current=nullptr;
};

template<typename TImage>
struct ImageRegionConstIteratorWithIndex : ImageRegionConstIterator<TImage> {
    typedef ImageRegionConstIterator<TImage> Superclass;
    typedef typename Superclass::PixelType PixelType;
    typedef typename Superclass::ImageType ImageType;
    typedef typename Superclass::ImageRegionType ImageRegionType;
    typedef typename ImageType::IndexType IndexType;
    typedef typename Superclass::ImagePointer ImagePointer;
    
    ImageRegionConstIteratorWithIndex(const ImageType * ptr,const ImageRegionType & region):Superclass(ptr,region){};
    
    void SetIndex(IndexType index){this->current = this->m_Ptr->GetBufferPointer()+this->m_Ptr->LinearIndex(index);}
    const IndexType GetIndex() const {
        IndexType result;
        ImageRegionType region = this->m_Ptr->GetBufferedRegion(); // ! Assumed to equal input region.
        const auto size = region.GetSize();
        const auto rIndex = region.GetIndex();
        DiscreteType linearIndex=this->current-this->m_Ptr->GetBufferPointer();
        for(int i=0; i<ImageType::ImageDimension; ++i){
            result[i]=rIndex[i]+(linearIndex%size[i]);
            linearIndex/=size[i];
        }
        return result;
    }
    std::string GetNameOfClass() const {return "ImageRegionConstIteratorWithIndex";}
};
    
    
template<typename TImage>
struct ImageRegionIteratorWithIndex : ImageRegionConstIteratorWithIndex<TImage> {
    typedef ImageRegionConstIteratorWithIndex<TImage> Superclass;
    typedef typename Superclass::ImageType ImageType;
    typedef typename Superclass::PixelType PixelType;
    typedef typename Superclass::ImagePointer ImagePointer;
    typedef typename Superclass::ImageRegionType ImageRegionType;
    ImageRegionIteratorWithIndex(ImageType * ptr,const ImageRegionType & region):Superclass(ptr,region){};
    PixelType & Value()       {return *const_cast<PixelType *>(this->current);} // Not for const iterator
};
}
#endif
