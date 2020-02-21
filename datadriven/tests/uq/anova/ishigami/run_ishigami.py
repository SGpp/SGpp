# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org
from argparse import ArgumentParser
from multiprocessing.process import Process
from uq.anova.ishigami.test_ishigami import run_ishigami_pce, \
    run_ishigami_sg, checkSobolIndices

if __name__ == '__main__':
    parallel = False

    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--table', default=False, action='store_true', help='just run the ones selected to come into the thesis')
    parser.add_argument('--out', default=False, action='store_true', help='write results to file')
    parser.add_argument('--plot', default=False, action='store_true', help='plot results (1d)')
    parser.add_argument('--surrogate', default="both", type=str, help="define which surrogate model should be used (sg, pce)")
    parser.add_argument('--parallel', default=False, action='store_true', help='run in parallel')
    args = parser.parse_args()

    scenarions_sg = {'gridType': []}
    scenarions_pce = {'sampler': []}
    if args.table:
        scenarions_sg = {'gridType': ["modPolyClenshawCurtis"],
                         'level': [2, 3, 4, 5],
                         'refinement': [None, 'var'],
                         'maxGridPoints': [65, 250]}

        scenarions_pce = {'sampler': ['fekete', 'leja'],
                         'degree': [5, 9]}
    else:
        scenarions_sg = {'gridType': ["linear",
                                      "poly",
                                      "modlinear",
                                      "modPolyClenshawCurtis",
                                      "modpoly"],
                         'level': [2, 3, 4, 5],
                         'refinement': ['var', 'simple',
                                        'exp', 'squared']}

        scenarions_pce = {'sampler': ['fekete', 'leja'],
                         'degree': [3, 5, 7, 9, 11]}

    processes = []
    if args.surrogate in ["both", "sg"]:
        for gridType in scenarions_sg['gridType']:
            for refinement in scenarions_sg['refinement']:
                maxGridPoints = [scenarions_sg['maxGridPoints'][-1]]
                if refinement is not None:
                    maxGridPoints = scenarions_sg['maxGridPoints']
                for level in scenarions_sg['level']:
                    for maxGridPoint in maxGridPoints:
                        print("-" * 80)
                        print("scenario: (%s, %i, %i, %s)" % (gridType, level, maxGridPoint, refinement))
                        print("-" * 80)
                        if args.parallel:
                            myargs = (gridType, level, maxGridPoint, False, refinement, args.out)
                            processes.append(Process(target=run_ishigami_sg, args=myargs))
                        else:
                            sobol_indices_analytic, sobol_indices, N = \
                                run_ishigami_sg(gridType,
                                                level,
                                                maxGridPoint,
                                                False,
                                                refinement,
                                                args.out)
                            if args.plot:
                                checkSobolIndices(sobol_indices_analytic, sobol_indices, N, args.plot)

    if args.surrogate in ["both", "pce"]:
        for sampler in scenarions_pce['sampler']:
            for degree_1d in scenarions_pce['degree']:
                print("-" * 80)
                print("scenario: (%s, %i)" % (sampler, degree_1d))
                print("-" * 80)
                if args.parallel:
                    myargs = (sampler, degree_1d, args.out)
                    processes.append(Process(target=run_ishigami_pce, args=myargs))
                else:
                    sobol_indices_analytic, sobol_indices, N = \
                        run_ishigami_pce(sampler, degree_1d, args.out)
                    if args.plot:
                        checkSobolIndices(sobol_indices_analytic, sobol_indices, N, args.plot)

    # run applications in parallel if there are any available
    for process in processes:
        process.start()
