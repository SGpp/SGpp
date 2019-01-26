import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from Grid import *
import time
from StandardCombi import *
from DimAdaptiveCombi import *


def performTestStandard(f, a, b, grid, lmin, maxLmax, dim, reference_solution):
    # calculate standard combination scheme results
    errorArrayStandard = []
    pointArray = []
    distinctFEvalArray = []
    standardCombi = StandardCombi(a, b, grid)
    for i in range(lmin + 1, lmin + maxLmax):
        scheme, error, result = standardCombi.perform_combi(lmin, i, f, reference_solution=reference_solution)
        errorArrayStandard.append(error / abs(reference_solution))
        pointArray.append(standardCombi.get_total_num_points())
        distinctFEvalArray.append(standardCombi.get_total_num_points(distinct_function_evals=True))
    return pointArray, distinctFEvalArray, errorArrayStandard


def performTestcaseArbitraryDim(f, a, b, adaptiveAlgorithmVector, maxtol, dim, maxLmax, grid=None, minLmin=1, maxLmin=3,
                                minTol=-1, doDimAdaptive=False):
    # realIntegral = scipy.integrate.dblquad(f, a, b, lambda x:a, lambda x:b, epsabs=1e-15, epsrel=1e-15)[0]
    reference_solution = f.getAnalyticSolutionIntegral(a, b)
    print("Exact integral", reference_solution)
    errorArray = []
    numEvaluationsArray = []
    numNaive = []
    numIdeal = []
    numFEvalIdeal = []
    # calculate refinements for different tolerance values
    for algorithm in adaptiveAlgorithmVector:
        errorArrayAlgorithm = []
        numEvaluationsArrayAlgorithm = []
        numNaiveAlgorithm = []
        numIdealAlgorithm = []
        numFEvalIdealAlgorithm = []
        '''
        for i in range(minTol, maxtol+1):
            start = time.time()
            if i == minTol:
                coarsening, combischeme, lmax, integral, numberOfEvaluations, error_array_new, num_point_array_new = algorithm[0].performSpatiallyAdaptiv(
                    algorithm[1], algorithm[2], f, algorithm[3], 10 ** -i, reference_solution=reference_solution)
                # errorArrayAlgorithm.append(abs(integral - realIntegral) / abs(realIntegral))
                errorArrayAlgorithm.extend(error_array_new)
                numEvaluationsArrayAlgorithm.append(numberOfEvaluations)
                # numIdealAlgorithm.extend(algorithm[0].get_total_num_points_arbitrary_dim(False))
                # numNaiveAlgorithm.append(algorithm[0].get_total_num_points_arbitrary_dim(True))
                numFEvalIdealAlgorithm.extend(num_point_array_new)
            else:
                if abs(integral - reference_solution) / abs(reference_solution) > 10 ** -i:
                    coarsening, combischeme, lmax, integral, numberOfEvaluations, error_array_new, num_point_array_new = algorithm[0].performSpatiallyAdaptiv(
                        algorithm[1], algorithm[2], f, algorithm[3], 10 ** -i, coarsening, reference_solution=reference_solution)
                    #errorArrayAlgorithm.append(abs(integral - realIntegral) / abs(realIntegral))
                    errorArrayAlgorithm.extend(error_array_new)
                    numEvaluationsArrayAlgorithm.append(numberOfEvaluations)
                    #numIdealAlgorithm.extend(algorithm[0].get_total_num_points_arbitrary_dim(False))
                    #numNaiveAlgorithm.append(algorithm[0].get_total_num_points_arbitrary_dim(True))
                    numFEvalIdealAlgorithm.extend(num_point_array_new)
            end = time.time()
            print("time spent in case", i, end - start)
            '''
        coarsening, combischeme, lmax, integral, numberOfEvaluations, error_array_new, num_point_array_new = algorithm[
            0].performSpatiallyAdaptiv(
            algorithm[1], algorithm[2], f, algorithm[3], 10 ** -maxtol, reference_solution=reference_solution)
        # errorArrayAlgorithm.append(abs(integral - realIntegral) / abs(realIntegral))
        errorArrayAlgorithm.extend(error_array_new)
        numEvaluationsArrayAlgorithm.append(numberOfEvaluations)
        # numIdealAlgorithm.extend(algorithm[0].get_total_num_points_arbitrary_dim(False))
        # numNaiveAlgorithm.append(algorithm[0].get_total_num_points_arbitrary_dim(True))
        numFEvalIdealAlgorithm.extend(num_point_array_new)

        errorArray.append(errorArrayAlgorithm)
        numEvaluationsArray.append(numEvaluationsArrayAlgorithm)
        numNaive.append(numNaiveAlgorithm)
        numIdeal.append(numIdealAlgorithm)
        numFEvalIdeal.append(numFEvalIdealAlgorithm)

    if doDimAdaptive:
        dimAdaptiveCombi = DimAdaptiveCombi(a,b, grid)
        scheme, error, result, errorArrayDimAdaptive, numFEvalIdealDimAdaptive = dimAdaptiveCombi.perform_combi(1,2,f,10**-maxtol, reference_solution=reference_solution)

    # calculate different standard combination scheme results
    xArrayStandard = []
    xFEvalArrayStandard = []
    errorArrayStandard = []
    for i in range(minLmin, maxLmin):
        xArrayStandardTest, xFEvalArrayStandardTest, errorArrayStandardTest = performTestStandard(f, a, b, grid, i,
                                                                                                  maxLmax - (i - 1) * (
                                                                                                              dim - 1),
                                                                                                  dim, reference_solution)
        xArrayStandard.append(xArrayStandardTest)
        xFEvalArrayStandard.append(xFEvalArrayStandardTest)
        errorArrayStandard.append(errorArrayStandardTest)
    # plot
    for i in range(maxLmin - minLmin):
        print(xArrayStandard[i], errorArrayStandard[i], "Number of Points Standard lmin= " + str(i + minLmin))
        print(xFEvalArrayStandard[i], errorArrayStandard[i], "Distinct f evaks Standard lmin= " + str(i + minLmin))
        # plt.loglog(xArrayStandard[i],errorArrayStandard[i],label='standardCombination lmin='+ str(i+minLmin))
        plt.loglog(xFEvalArrayStandard[i], errorArrayStandard[i],
                   label='standardCombination distinct f evals lmin=' + str(i + minLmin))
    if doDimAdaptive:    
        print(numFEvalIdealDimAdaptive, errorArrayDimAdaptive, "Number of Points DimAdaptive lmin= 1")
        plt.loglog(numFEvalIdealDimAdaptive, errorArrayDimAdaptive, label="Number of Points DimAdaptive lmin= 1")
    for i in range(len(adaptiveAlgorithmVector)):
        #print(numNaive[i], errorArray[i], adaptiveAlgorithmVector[i][4] + ' Naive evaluation')
        #print(numIdeal[i], errorArray[i], adaptiveAlgorithmVector[i][4] + ' total points')
        print(numFEvalIdeal[i], errorArray[i], adaptiveAlgorithmVector[i][4] + ' distinct f evals')
        # plt.loglog(numNaive[i],errorArray[i],label= adaptiveAlgorithmVector[i][3] +' Naive evaluation')
        # plt.loglog(numIdeal[i],errorArray[i],label=adaptiveAlgorithmVector[i][3] +' total points')
        plt.loglog(numFEvalIdeal[i], errorArray[i], label=adaptiveAlgorithmVector[i][4] + ' distinct f evals')
    plt.legend(bbox_to_anchor=(3, 1), loc="upper right")
    plt.xlabel('Number of points')
    plt.ylabel('Approximation error')
    # plt.savefig('convergence.pdf', bbox_inches='tight')
    plt.show()
