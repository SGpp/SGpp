Cube Data Set
=========================
From:
Dirk Pflueger (pflueged@in.tum.de)

Description:
A d-dimensional dataset, with uniformly distributed points (random) in
[0,1]^d. 
All points with class -1 are in a subcube with volume 1/2, originating
at the origin.

Files:
create_dataset.py
  Creates datasets as whitespace separated list
create_scripts.py
  Creates script for converter.py to convert files into arff-format

The scripts that created the data sets:
create_1000tr500te.sh
  Created by:
  > python create_scripts.py 1 > create_1000tr500te.sh
dat2arff.sh
  Created by:
  > python create_scripts.py 2 > dat2arff.sh

The data files are not included but have to be created manually, as
they take too much storage space.
Update: Now in .gz-Format!
