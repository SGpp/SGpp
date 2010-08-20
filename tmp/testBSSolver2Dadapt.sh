#!/bin/sh

LEVELS="6 7 8 9 10"
REFINEMENTS="1 2 3 4 5 6 7"
MAXLEVELS="10 11 12 13 14"
REFINETHRESHOLDS="1e-1 5e-2 1e-2 5e-3 1e-3 5e-4 1e-4 5e-5 1e-5 5e-6"
ADAPTSOLVEMODE="coarsen coarsenNrefine"
COARSENTHRESHOLD="1e-6 5e-7 1e-7 5e-8 1e-8 5e-9 1e-9"

#Do tests without adaptivity during solving with classic adativity
for level in $LEVELS
do
	for refNums in $REFINEMENTS
	do
		for refThres in $REFINETHRESHOLDS
		do
			OUTFILE="BS_2d_log_put_adapt_classic_l-${level}_nr-${refNums}_tr-${refThres}.out"
			echo "./BSSolver solveNDadaptSurplus log 2 $level BStest2DMC.log.bound BStest2DMC.stoch 1.0 std_euro_put 0.00 1.0 0.05 CrNic 10000 0.00001 classic 0 $refNums $refThres none 1e-6"
			./BSSolver solveNDadaptSurplus log 2 $level BStest2DMC.log.bound BStest2DMC.stoch 1.0 std_euro_put 0.00 1.0 0.05 CrNic 10000 0.00001 classic 0 $refNums $refThres none 1e-6 > $OUTFILE
		done
	done
done

#Do tests without adaptivity during solving with maxLevel adativity
for level in $LEVELS
do
	for refNums in $REFINEMENTS
	do
		for refThres in $REFINETHRESHOLDS
		do
			for maxlevel in $MAXLEVELS
			do
				OUTFILE="BS_2d_log_put_adapt_maxLevel_l-${level}_nr-${refNums}_tr-${refThres}_mlr-${maxlevel}.out"
				echo "./BSSolver solveNDadaptSurplus log 2 $level BStest2DMC.log.bound BStest2DMC.stoch 1.0 std_euro_put 0.00 1.0 0.05 CrNic 10000 0.00001 maxLevel $maxlevel $refNums $refThres none 1e-6"
				./BSSolver solveNDadaptSurplus log 2 $level BStest2DMC.log.bound BStest2DMC.stoch 1.0 std_euro_put 0.00 1.0 0.05 CrNic 10000 0.00001 maxLevel $maxlevel $refNums $refThres none 1e-6 > $OUTFILE
			done
		done
	done
done

#Do tests with adaptivity during solving with classic adativity
for level in $LEVELS
do
	for refNums in $REFINEMENTS
	do
		for refThres in $REFINETHRESHOLDS
		do
			for adaptsolve in $ADAPTSOLVEMODE
			do
				for coaThres in $COARSENTHRESHOLD
				do
					OUTFILE="BS_2d_log_put_adapt_classic_l-${level}_nr-${refNums}_tr-${refThres}_as-${adaptsolve}_tc-${coaThres}.out"
					echo "./BSSolver solveNDadaptSurplus log 2 $level BStest2DMC.log.bound BStest2DMC.stoch 1.0 std_euro_put 0.00 1.0 0.05 CrNic 10000 0.00001 classic 0 $refNums $refThres $adaptsolve $coaThres"
					./BSSolver solveNDadaptSurplus log 2 $level BStest2DMC.log.bound BStest2DMC.stoch 1.0 std_euro_put 0.00 1.0 0.05 CrNic 10000 0.00001 classic 0 $refNums $refThres $adaptsolve $coaThres > $OUTFILE
				done
			done
		done
	done
done

#Do tests with adaptivity during solving with maxLevel adativity
for level in $LEVELS
do
	for refNums in $REFINEMENTS
	do
		for refThres in $REFINETHRESHOLDS
		do
			for adaptsolve in $ADAPTSOLVEMODE
			do
				for coaThres in $COARSENTHRESHOLD
				do
					for maxlevel in $MAXLEVELS
					do
						OUTFILE="BS_2d_log_put_adapt_maxLevel_l-${level}_nr-${refNums}_tr-${refThres}_mlr-${maxlevel}_as-${adaptsolve}_tc-${coaThres}.out"
						echo "./BSSolver solveNDadaptSurplus log 2 $level BStest2DMC.log.bound BStest2DMC.stoch 1.0 std_euro_put 0.00 1.0 0.05 CrNic 10000 0.00001 maxLevel $maxlevel $refNums $refThres $adaptsolve $coaThres"
						./BSSolver solveNDadaptSurplus log 2 $level BStest2DMC.log.bound BStest2DMC.stoch 1.0 std_euro_put 0.00 1.0 0.05 CrNic 10000 0.00001 maxLevel $maxlevel $refNums $refThres $adaptsolve $coaThres > $OUTFILE
					done
				done
			done
		done
	done
done