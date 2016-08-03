all:
	nvcc -arch=sm_35 -dc kernels.cu -rdc=true
	nvcc -arch=sm_35 -dlink -o kernels_l.o kernels.o -lcudart
	nvcc -arch=sm_35 -dc MortonOrderKernel.cu -rdc=true
	nvcc -arch=sm_35 -dlink -o MortonOrderKernel_l.o MortonOrderKernel.o -lcudart
	nvcc -arch=sm_35 -dc MultiEvalKernel.cu -rdc=true
	nvcc -arch=sm_35 -dlink -o MultiEvalKernel_l.o MultiEvalKernel.o -lcudart
	g++ -Wall -c -std=c++11 MortonOrder.cpp -I /home/malte/Documents/Uni/HiWi/SGppCD16/SGpp/base/src/ -L /home/malte/Documents/Uni/HiWi/SGppCD16/SGpp/lib/sgpp/ -I /usr/local/cuda/include -L /usr/local/cuda/lib64/ -lsgppbase -lcudart
	g++ -Wall -c -std=c++11 OperationMultiEvalCuda.cpp -I /home/malte/Documents/Uni/HiWi/SGppCD16/SGpp/base/src/ -L /home/malte/Documents/Uni/HiWi/SGppCD16/SGpp/lib/sgpp/ -I /usr/local/cuda/include -L /usr/local/cuda/lib64/ -lsgppbase -lcudart  
	ar rvs libCudaOpEval.a *.o
	
test: test.cpp libCudaOpEval.a
	rm -f test
	g++ -Wall -o test -std=c++11 test.cpp -I /home/malte/Documents/Uni/HiWi/SGppCD16/SGpp/base/src/ -I /home/malte/Documents/Uni/HiWi/SGppCD16/SGpp/datadriven/src -L /home/malte/Documents/Uni/HiWi/SGppCD16/SGpp/lib/sgpp/ -I /usr/local/cuda/include -L /usr/local/cuda/lib64/ -L ./ -lCudaOpEval -lcudart -lcuda -lsgppbase -lsgppdatadriven -lsgppsolver -lsgpppde

clean:
	rm *.o
	rm *.a
	rm test
