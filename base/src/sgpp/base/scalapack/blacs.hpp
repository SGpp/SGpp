#pragma once

extern "C" {

void blacs_gridinfo_(int &icontxt, int &nprow, int &npcol, int &myprow, int &mypcol);

void blacs_gridexit_(int &icontxt);

void blacs_exit_(int &cont);
}