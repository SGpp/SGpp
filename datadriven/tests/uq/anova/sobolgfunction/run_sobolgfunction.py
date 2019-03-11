'''
Created on Sep 19, 2016

@author: franzefn
'''
from argparse import ArgumentParser
from multiprocessing.process import Process
from uq.anova.sobolgfunction.test_sobolgfunction import run_sobol_g_function_pce, \
    run_sobol_g_function_sg, checkSobolIndices

if __name__ == '__main__':
    parallel = False

    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--table', default=False, action='store_true', help='just run the ones selected to come into the thesis')
    parser.add_argument('--surrogate', default="both", type=str, help="define which surrogate model should be used (sg, pce)")
    parser.add_argument('--full', default=False, action='store_true', help='run the full model')
    parser.add_argument('--reduced', default=False, action='store_true', help='run the reduced model')
    parser.add_argument('--out', default=False, action='store_true', help='write results to file')
    parser.add_argument('--plot', default=False, action='store_true', help='plot results (1d)')
    parser.add_argument('--parallel', default=False, action='store_true', help='run in parallel')
    args = parser.parse_args()

    scenarions_sg = {'gridType': []}
    scenarions_pce = {'sampler': []}
    if args.table:
        scenarions_sg = {'full': {'gridType': ["modPolyClenshawCurtis"],
                                  'level': [2, 3],
                                  'refinement': [None, 'var'],
                                  'maxGridPoints': [70, 140]},
                         'reduced': {'gridType': ["modPolyClenshawCurtis"],
                                     'level': [2, 3, 4],
                                     'refinement': [None, 'var'],
                                     'maxGridPoints': [100, 200]},
                                     }

        scenarions_pce = {'full': {'sampler': ['fekete', 'leja'],
                                   'degree': [2]},
                          'reduced': {'sampler': ['fekete', 'leja'],
                                   'degree': [3, 5, 7]}}
    else:
        scenarions_sg = {'full': {'gridType': ["linear",
                                               "poly",
                                               "modlinear",
                                               "modPolyClenshawCurtis",
                                               "modpoly"],
                                  'level': [2, 3],
                                  'refinement': ['var', 'simple',
                                                 'exp', 'squared']},
                         'reduced': {'gridType': ["linear",
                                                  "poly",
                                                  "modlinear",
                                                  "modPolyClenshawCurtis",
                                                  "modpoly"],
                                     'level': [2],
                                     'refinement': ['var', 'simple',
                                                    'exp', 'squared']}}

        scenarions_pce = {'full': {'sampler': ['fekete', 'leja'],
                                  'degree': [2, 3, 4]},
                         'reduced': {'sampler': ['fekete', 'leja'],
                                     'degree': [2, 3, 4, 5, 6, 7, 8]}}

    processes = []
    if args.surrogate in ["both", "sg"]:
        for model, surrogate in list(scenarions_sg.items()):
            if args.full and model == "reduced":
                continue
            if args.reduced and model == "full":
                continue

            for gridType in surrogate['gridType']:
                for refinement in surrogate['refinement']:
                    maxGridPoints = [surrogate['maxGridPoints'][-1]]
                    if refinement is not None:
                        maxGridPoints = surrogate['maxGridPoints']
                    for level in surrogate['level']:
                        for maxGridPoint in maxGridPoints:
                            print("-" * 80)
                            print("scenario: (%s, %s, %i, %i, %s)" % (model, gridType, level, maxGridPoint, refinement))
                            print("-" * 80)

                            if args.parallel:
                                myargs = (model == "full",
                                          gridType, level, maxGridPoint, False, refinement, args.out)
                                processes.append(Process(target=run_sobol_g_function_sg, args=myargs))
                            else:
                                sobol_indices_analytic, sobol_indices, N = \
                                    run_sobol_g_function_sg(model == "full",
                                                            gridType,
                                                            level,
                                                            maxGridPoint,
                                                            False,
                                                            refinement,
                                                            args.out)
                                if args.plot:
                                    checkSobolIndices(sobol_indices_analytic, sobol_indices, N, False)

    if args.surrogate in ["both", "pce"]:
        for model, surrogate in list(scenarions_pce.items()):
            if args.full and model == "reduced":
                continue
            if args.reduced and model == "full":
                continue

            for sampler in surrogate['sampler']:
                for degree_1d in surrogate['degree']:
                    print("-" * 80)
                    print("scenario: (%s, %s, %i)" % (model, sampler, degree_1d))
                    print("-" * 80)

                    if args.parallel:
                        myargs = (model == "full", sampler, degree_1d, True)
                        processes.append(Process(target=run_sobol_g_function_pce, args=myargs))
                    else:
                        sobol_indices_analytic, sobol_indices, N = \
                            run_sobol_g_function_pce(model == "full", sampler, degree_1d, True)

                        if args.plot:
                            checkSobolIndices(sobol_indices_analytic, sobol_indices, N, False)

    # run applications in parallel if there are any available
    for process in processes:
        process.start()
