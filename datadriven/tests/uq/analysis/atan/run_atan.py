'''
Created on Sep 19, 2016

@author: franzefn
'''
from argparse import ArgumentParser
from multiprocessing.process import Process
from uq.analysis.atan.test_atan import run_atan_pce, run_atan_sg

if __name__ == '__main__':
    parallel = False

    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--surrogate', default="both", type=str, help="define which surrogate model should be used (sg, pce)")
    parser.add_argument('--model', default="beta", type=str, help="define which probabilistic model should be used")
    parser.add_argument('--out', default=False, action='store_true', help='write results to file')
    parser.add_argument('--plot', default=False, action='store_true', help='plot functions (2d)')
    parser.add_argument('--parallel', default=False, action='store_true', help='run in parallel')
    args = parser.parse_args()

    scenarions_sg = {'uniform': {'gridType': ["linearBoundary",
                                              "modlinear",
                                              "linearClenshawCurtisBoundary",
                                              "modLinearClenshawCurtis",
                                              "polyBoundary",
                                              "modpoly",
                                              "polyClenshawCurtisBoundary",
                                              "modPolyClenshawCurtis"],
                                 'full': [True, False],
                                 'level': [2, 3, 4],
                                 'refinement': [
                                                None,
                                                'simple',
                                                'weighted',
                                                'squared',
                                                'mean_squared',
                                                'var',
                                                'exp',
                                                'l2'],
                                 'maxGridPoints': [3000]},
                     'beta': {'gridType': [
                                           "linearBoundary",
#                                            "modlinear",
#                                            "linearClenshawCurtisBoundary",
#                                            "modLinearClenshawCurtis",
                                           "polyBoundary",
#                                            "modpoly",
#                                            "polyClenshawCurtisBoundary",
#                                            "modPolyClenshawCurtis"
                                           ],
                               'full': [True, False],
                               'level': [2],  # , 3, 4],
                               'refinement': [
                                              None,
                                              'simple',
                                              'weighted',
                                              'mean_squared',
                                              'squared',
                                              'var',
                                              'exp',
                                              'l2'],
                               'maxGridPoints': [3000]}}

    scenarions_pce = {'uniform': {'sampler': ['gauss', 'gauss_leja', 'fekete', 'leja', ],
                                  'expansion': ["full_tensor", "total_degree"],
                                  'max_num_samples': [4000]},
                      'beta': {'sampler': ['gauss', 'gauss_leja', 'fekete', 'leja', ],
                               'expansion': ["full_tensor", "total_degree"],
                               'max_num_samples': [4000]}}

    processes = []
    if args.surrogate in ["both", "sg"]:
        surrogate = scenarions_sg[args.model]
        for level in surrogate["level"]:
            for gridType in surrogate['gridType']:
                for refinement in surrogate['refinement']:
                    maxGridPoints = surrogate['maxGridPoints']
                    fulls = surrogate['full']
                    if refinement is not None:
                        fulls = [False]
                        maxGridPoints = [surrogate['maxGridPoints'][-1]]

                    for full in fulls:
                        for maxGridPoint in maxGridPoints:
                            for full in fulls:
                                print("-" * 80)
                                print("scenario: (%s, %s, %i, %s, %s)" % (args.model, gridType, maxGridPoint, refinement,
                                                                          "full" if full else "sparse"))
                                print("-" * 80)

                                if args.parallel:
                                    myargs = (args.model, gridType, level, maxGridPoint,
                                              1, full, refinement, args.out,
                                              args.plot)
                                    processes.append(Process(target=run_atan_sg, args=myargs))
                                else:
                                    run_atan_sg(args.model,
                                                gridType,
                                                level,
                                                maxGridPoint,
                                                1,
                                                full,
                                                refinement,
                                                args.out,
                                                args.plot)

    if args.surrogate in ["both", "pce"]:
        surrogate = scenarions_pce[args.model]
        for sampler in surrogate['sampler']:
            for expansion in surrogate["expansion"]:
                for max_num_samples in surrogate["max_num_samples"]:
                    print("-" * 80)
                    print("scenario: (%s, %s, %s, %i)" % (args.model, sampler, expansion, max_num_samples))
                    print("-" * 80)

                    if args.parallel:
                        myargs = (args.model, sampler, expansion, max_num_samples, args.out, args.plot)
                        processes.append(Process(target=run_atan_pce, args=myargs))
                    else:
                        run_atan_pce(args.model, sampler, expansion, max_num_samples, args.out, args.plot)

    # run applications in parallel if there are any available
    for process in processes:
        process.start()
