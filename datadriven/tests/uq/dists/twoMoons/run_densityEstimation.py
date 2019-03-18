'''
Created on Sep 19, 2016

@author: franzefn
'''
from argparse import ArgumentParser
from multiprocessing.process import Process
from estimateDensity import density_configs, run_densityEstimation

if __name__ == '__main__':
    parallel = False

    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--out', default=False, action='store_true', help='write results to file')
    parser.add_argument('--plot', default=False, action='store_true', help='do some plotting while estimation')
    parser.add_argument('--parallel', default=False, action='store_true', help='run in parallel')
    args = parser.parse_args()

    scenarions = {'density_configs': ["sgde_zero"], #density_configs,
                  "function_configs": [("two_moons", 2), ("mult_beta", 2)]}
    processes = []
    for estimationMethod in scenarions['density_configs']:
        for functionName, numDims in scenarions["function_configs"]:
            print("-" * 80)
            print("scenario: %s, d=%i, %s" % (functionName, numDims, estimationMethod))
            print("-" * 80)

            if args.parallel:
                myargs = (functionName,
                          estimationMethod,
                          10,
                          numDims,
                          20000,
                          "join",
                          "ml",
                          args.out,
                          args.plot)
                processes.append(Process(target=run_densityEstimation, args=myargs))
            else:
                run_densityEstimation(functionName,
                                      estimationMethod,
                                      kfold=10,
                                      numDims=numDims,
                                      numSamples=20000,
                                      candidates="join",
                                      bandwidthOptimizationType="ml",
                                      out=args.out,
                                      plot=args.plot)

    # run applications in parallel if there are any available
    for process in processes:
        process.start()
