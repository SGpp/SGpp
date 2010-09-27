comp_comand="g++ -O3 -I../src/sgpp -L../lib/sgpp -o $1 $1.cpp -lsgpp"
echo $comp_comand
$comp_comand