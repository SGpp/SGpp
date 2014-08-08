#!/bin/sed -f 

s/@attribute/@ATTRIBUTE/
s/@relation/@RELATION/
s/@data/@DATA/
s/'-1'$/-1/
s/'1'$/1/
