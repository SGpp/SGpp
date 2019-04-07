Friedman1-3 Data Set
=========================
From:
Friedman, J. H. Multivariate adaptive regression splines. Ann. Statis., 19, 1-131. 1991.

friedman.py generates Friedman1-3 datasets, e.g.

# create files
# note: to generate comparable files, use seed 123456 for training, 234567 for testing, and 345678 for validation
# note: csv-files are created if the output file name is specified and contains ".csv", otherwise the output format is arff

python friedman.py -f 1 -N 10000 -o friedman1_10000_train.arff.gz [--seed 123456]
python friedman.py -f 1 -N 10000 -o friedman1_10000_test.arff.gz  [--seed 234567]
python friedman.py -f 1 -N 90000 -o friedman1_90000_train.arff.gz [--seed 123456]
								  
python friedman.py -f 2 -N 10000 -o friedman2_10000_train.arff.gz [--seed 123456]
python friedman.py -f 2 -N 10000 -o friedman2_10000_test.arff.gz  [--seed 234567]
python friedman.py -f 2 -N 90000 -o friedman2_90000_train.arff.gz [--seed 123456]
								  
python friedman.py -f 3 -N 10000 -o friedman3_10000_train.arff.gz [--seed 123456]
python friedman.py -f 3 -N 10000 -o friedman3_10000_test.arff.gz  [--seed 234567]
python friedman.py -f 3 -N 90000 -o friedman3_90000_train.arff.gz [--seed 123456]


# normalization
python ../../bin/converter.py --noclassnormalization --min 0 --max 100 --min 125.66 --max 1859.3 --min 0 --max 1 --min 1 --max 11 -i friedman2_10000_train.arff.gz -i friedman2_10000_test.arff.gz -i friedman2_90000_train.arff.gz -o friedman2_10000_train_norm.arff.gz -o friedman2_10000_test_norm.arff.gz -o friedman2_90000_train_norm.arff.gz 

python ../../bin/converter.py --noclassnormalization --min 0 --max 100 --min 125.66 --max 1859.3 --min 0 --max 1 --min 1 --max 11 -i friedman3_10000_train.arff.gz -i friedman3_10000_test.arff.gz -i friedman3_90000_train.arff.gz -o friedman3_10000_train_norm.arff.gz -o friedman3_10000_test_norm.arff.gz -o friedman3_90000_train_norm.arff.gz 

