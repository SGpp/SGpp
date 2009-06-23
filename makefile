CC = g++
DEBUG = 
SRCDIR = src/sgpp/
CFLAGS = -Wall -pedantic -ansi -c -Wno-long-long -fno-strict-aliasing -fopenmp -O3 -funroll-loops -pthread -DUSEOMP -Isrc/sgpp/  -Wno-deprecated
LFLAGS = -Wall -pedantic -ansi -fopenmp -O3 -pthread -Wno-deprecated
OBJS = DataVector.o Grid.o LinearBoundaryGrid.o LinearBoundaryUScaledGrid.o LinearGrid.o ModLinearGrid.o ModPolyGrid.o PolyGrid.o StandardGridGenerator.o BoundaryGridGenerator.o BoundaryUScaledGridGenerator.o OperationBLinear.o OperationEvalLinear.o OperationHierarchisationLinear.o OperationBModLinear.o OperationEvalModLinear.o OperationHierarchisationModLinear.o OperationBModPoly.o OperationEvalModPoly.o OperationHierarchisationModPoly.o OperationBPoly.o OperationEvalPoly.o OperationHierarchisationPoly.o OperationHierarchisationLinearBoundaryUScaled.o OperationEvalLinearBoundaryUScaled.o OperationBLinearBoundaryUScaled.o OperationHierarchisationLinearBoundary.o OperationEvalLinearBoundary.o OperationBLinearBoundary.o ARFFTools.o ApplyDMMatrix.o Classifier.o NativeCppClassifier.o
OBJSSRC = tmp/build_native/DataVector.o tmp/build_native/Grid.o tmp/build_native/LinearBoundaryGrid.o tmp/build_native/LinearBoundaryUScaledGrid.o tmp/build_native/LinearGrid.o tmp/build_native/ModLinearGrid.o tmp/build_native/ModPolyGrid.o tmp/build_native/PolyGrid.o tmp/build_native/StandardGridGenerator.o tmp/build_native/BoundaryGridGenerator.o tmp/build_native/BoundaryUScaledGridGenerator.o tmp/build_native/OperationBLinear.o tmp/build_native/OperationEvalLinear.o tmp/build_native/OperationHierarchisationLinear.o tmp/build_native/OperationBModLinear.o tmp/build_native/OperationEvalModLinear.o tmp/build_native/OperationHierarchisationModLinear.o tmp/build_native/OperationBModPoly.o tmp/build_native/OperationEvalModPoly.o tmp/build_native/OperationHierarchisationModPoly.o tmp/build_native/OperationBPoly.o tmp/build_native/OperationEvalPoly.o tmp/build_native/OperationHierarchisationPoly.o tmp/build_native/OperationHierarchisationLinearBoundaryUScaled.o tmp/build_native/OperationEvalLinearBoundaryUScaled.o tmp/build_native/OperationBLinearBoundaryUScaled.o tmp/build_native/OperationHierarchisationLinearBoundary.o tmp/build_native/OperationEvalLinearBoundary.o tmp/build_native/OperationBLinearBoundary.o tmp/build_native/ARFFTools.o tmp/build_native/ApplyDMMatrix.o tmp/build_native/Classifier.o tmp/build_native/NativeCppClassifier.o

NativeCppClassifier : $(OBJS)
	mkdir -p tmp/build_native
	\mv *.o tmp/build_native/
	$(CC) $(LFLAGS) $(OBJSSRC) -o bin/NativeCppClassifier

DataVector.o : src/sgpp/data/DataVector.cpp
	$(CC) $(CFLAGS) src/sgpp/data/DataVector.cpp

Grid.o : src/sgpp/grid/Grid.cpp
	$(CC) $(CFLAGS) src/sgpp/grid/Grid.cpp

LinearBoundaryGrid.o : src/sgpp/grid/type/LinearBoundaryGrid.cpp
	$(CC) $(CFLAGS) src/sgpp/grid/type/LinearBoundaryGrid.cpp

LinearBoundaryUScaledGrid.o : src/sgpp/grid/type/LinearBoundaryUScaledGrid.cpp
	$(CC) $(CFLAGS) src/sgpp/grid/type/LinearBoundaryUScaledGrid.cpp
	
LinearGrid.o : src/sgpp/grid/type/LinearGrid.cpp
	$(CC) $(CFLAGS) src/sgpp/grid/type/LinearGrid.cpp
	
ModLinearGrid.o : src/sgpp/grid/type/ModLinearGrid.cpp
	$(CC) $(CFLAGS) src/sgpp/grid/type/ModLinearGrid.cpp

ModPolyGrid.o : src/sgpp/grid/type/ModPolyGrid.cpp
	$(CC) $(CFLAGS) src/sgpp/grid/type/ModPolyGrid.cpp

PolyGrid.o : src/sgpp/grid/type/PolyGrid.cpp
	$(CC) $(CFLAGS) src/sgpp/grid/type/PolyGrid.cpp
	
StandardGridGenerator.o : src/sgpp/grid/generation/StandardGridGenerator.cpp
	$(CC) $(CFLAGS) src/sgpp/grid/generation/StandardGridGenerator.cpp
	
BoundaryGridGenerator.o : src/sgpp/grid/generation/BoundaryGridGenerator.cpp
	$(CC) $(CFLAGS) src/sgpp/grid/generation/BoundaryGridGenerator.cpp

BoundaryUScaledGridGenerator.o : src/sgpp/grid/generation/BoundaryUScaledGridGenerator.cpp
	$(CC) $(CFLAGS) src/sgpp/grid/generation/BoundaryUScaledGridGenerator.cpp

OperationBLinear.o : src/sgpp/basis/linear/operation/OperationBLinear.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/linear/operation/OperationBLinear.cpp

OperationEvalLinear.o : src/sgpp/basis/linear/operation/OperationEvalLinear.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/linear/operation/OperationEvalLinear.cpp

OperationHierarchisationLinear.o : src/sgpp/basis/linear/operation/OperationHierarchisationLinear.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/linear/operation/OperationHierarchisationLinear.cpp

OperationBModLinear.o : src/sgpp/basis/modlinear/operation/OperationBModLinear.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/modlinear/operation/OperationBModLinear.cpp

OperationEvalModLinear.o : src/sgpp/basis/modlinear/operation/OperationEvalModLinear.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/modlinear/operation/OperationEvalModLinear.cpp

OperationHierarchisationModLinear.o : src/sgpp/basis/modlinear/operation/OperationHierarchisationModLinear.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/modlinear/operation/OperationHierarchisationModLinear.cpp

OperationBModPoly.o : src/sgpp/basis/modpoly/operation/OperationBModPoly.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/modpoly/operation/OperationBModPoly.cpp

OperationEvalModPoly.o : src/sgpp/basis/modpoly/operation/OperationEvalModPoly.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/modpoly/operation/OperationEvalModPoly.cpp

OperationHierarchisationModPoly.o : src/sgpp/basis/modpoly/operation/OperationHierarchisationModPoly.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/modpoly/operation/OperationHierarchisationModPoly.cpp

OperationBPoly.o : src/sgpp/basis/poly/operation/OperationBPoly.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/poly/operation/OperationBPoly.cpp

OperationEvalPoly.o : src/sgpp/basis/poly/operation/OperationEvalPoly.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/poly/operation/OperationEvalPoly.cpp

OperationHierarchisationPoly.o : src/sgpp/basis/poly/operation/OperationHierarchisationPoly.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/poly/operation/OperationHierarchisationPoly.cpp

OperationHierarchisationLinearBoundaryUScaled.o : src/sgpp/basis/linearboundaryUScaled/operation/OperationHierarchisationLinearBoundaryUScaled.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/linearboundaryUScaled/operation/OperationHierarchisationLinearBoundaryUScaled.cpp

OperationEvalLinearBoundaryUScaled.o : src/sgpp/basis/linearboundaryUScaled/operation/OperationEvalLinearBoundaryUScaled.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/linearboundaryUScaled/operation/OperationEvalLinearBoundaryUScaled.cpp

OperationBLinearBoundaryUScaled.o : src/sgpp/basis/linearboundaryUScaled/operation/OperationBLinearBoundaryUScaled.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/linearboundaryUScaled/operation/OperationBLinearBoundaryUScaled.cpp

OperationHierarchisationLinearBoundary.o : src/sgpp/basis/linearboundary/operation/OperationHierarchisationLinearBoundary.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/linearboundary/operation/OperationHierarchisationLinearBoundary.cpp

OperationEvalLinearBoundary.o : src/sgpp/basis/linearboundary/operation/OperationEvalLinearBoundary.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/linearboundary/operation/OperationEvalLinearBoundary.cpp

OperationBLinearBoundary.o : src/sgpp/basis/linearboundary/operation/OperationBLinearBoundary.cpp
	$(CC) $(CFLAGS) src/sgpp/basis/linearboundary/operation/OperationBLinearBoundary.cpp
	
ARFFTools.o : src/sgpp/tools/classification/ARFFTools.cpp
	$(CC) $(CFLAGS) src/sgpp/tools/classification/ARFFTools.cpp
	
ApplyDMMatrix.o : src/sgpp/solver/cg/ApplyMatrix/classification/ApplyDMMatrix.cpp
	$(CC) $(CFLAGS) src/sgpp/solver/cg/ApplyMatrix/classification/ApplyDMMatrix.cpp

Classifier.o : src/sgpp/application/classification/Classifier.cpp
	$(CC) $(CFLAGS) src/sgpp/application/classification/Classifier.cpp
	
NativeCppClassifier.o : src/sgpp/NativeCppClassifier.cpp
	$(CC) $(CFLAGS) src/sgpp/NativeCppClassifier.cpp
	
clean:
	\rm tmp/build_native/*.o bin/NativeCppClassifier *.o
