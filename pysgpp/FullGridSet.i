// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

%include "carrays.i"
%array_class(size_t, size_tArray);

%feature("notabstract") FullGridSet;
%rename(__len__) FullGridSet::getSize;
%rename(__getitem__) FullGridSet::get(unsigned int index,unsigned int pos);
%rename(__setitem__) FullGridSet::set(unsigned int index,unsigned int pos, double val );

namespace sg{
class FullGridSet{
    FullGridSet(size_t dimension,size_t n);
    FullGridSet(size_t dimension,size_t *n,const char *type);   
    FULLGRID* getGrids();
    FullGrid* at(size_t index);
    size_t getSize();   
    double get(size_t index,unsigned int pos);
    void set(size_t index,unsigned int pos, double val );
    size_t sizeOf(unsigned int index);
    void deCompose(GridStorage *storage,DataVector alpha);
    void reCompose(GridStorage *storage,DataVector *alpha);   
    FullGridSet(size_t dimension,size_t n,const char *type);
    BoundingBox* getBoundingBox();
    void getCoefs(DataVector &v);
};
}