comp_comand="g++ -O0 -g3 -D_GLIBCXX_DEBUG -I../src/sgpp -L../lib/sgpp -o $1 $1.cpp -lsgpp"
echo $comp_comand
$comp_comand