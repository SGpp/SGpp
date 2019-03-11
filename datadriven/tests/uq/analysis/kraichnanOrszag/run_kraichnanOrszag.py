'''
Created on Sep 19, 2016

@author: franzefn
'''
from argparse import ArgumentParser
from multiprocessing.process import Process
from uq.analysis.kraichnanOrszag.test_kraichnanOrszag import run_kraichnanOrszag_sg, run_kraichnanOrszag_mc

if __name__ == '__main__':
    parallel = False

    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--surrogate', default="both", type=str, help="define which surrogate model should be used (sg, mc)")
    parser.add_argument('--model', default="uniform", type=str, help="define which probabilistic model should be used")
    parser.add_argument('--out', default=False, action='store_true', help='write results to file')
    parser.add_argument('--parallel', default=False, action='store_true', help='run in parallel')
    args = parser.parse_args()

    scenarions_sg = {'uniform': {'gridType': ["linearBoundary",
                                              "modlinear",
                                              "polyBoundary",
                                              "modpoly"],
                                 'settings': [3],
                                 'qois': ["y1", "y3"],
                                 'level': [3],
                                 'refinement': [
                                                None,
                                                'simple',
                                                'exp',
                                                'l2'],
                                 'maxGridPoints': [1000]}}

    processes = []

    if args.surrogate in ["both", "sg"]:
        surrogate = scenarions_sg[args.model]
        maxGridPoints = surrogate['maxGridPoints'][0]
        for setting in surrogate["settings"]:
            for qoi in surrogate["qois"]:
                for level in surrogate["level"]:
                    for gridType in surrogate['gridType']:
                        for refinement in surrogate['refinement']:
                            print("-" * 80)
                            print("scenario: SG (%s, %s, %i, %s)" % (args.model,
                                                                     gridType,
                                                                     maxGridPoints,
                                                                     refinement))
                            print("-" * 80)

                            if args.parallel:
                                myargs = (gridType,
                                          level,
                                          maxGridPoints,
                                          1,
                                          False,
                                          refinement,
                                          setting,
                                          qoi,
                                          args.out,
                                          False)
                                processes.append(Process(target=run_kraichnanOrszag_sg,
                                                         args=myargs))
                            else:
                                run_kraichnanOrszag_sg(gridType,
                                                       level,
                                                       maxGridPoints,
                                                       1,
                                                       False,
                                                       refinement,
                                                       setting,
                                                       qoi,
                                                       args.out,
                                                       False)

    if args.surrogate in ["both", "mc"]:
        surrogate = scenarions_sg[args.model]
        maxSamples = 100000
        for setting in surrogate["settings"]:
            for qoi in surrogate["qois"]:
                print("-" * 80)
                print("scenario: MC (%s, %i)" % (args.model,
                                                 maxSamples))
                print("-" * 80)

                if args.parallel:
                    myargs = (setting,
                              qoi,
                              maxSamples,
                              args.out,
                              False)
                    processes.append(Process(target=run_kraichnanOrszag_mc,
                                             args=myargs))
                else:
                    run_kraichnanOrszag_mc(setting,
                                           qoi,
                                           maxSamples,
                                           args.out,
                                           False)


    # run applications in parallel if there are any available
    for process in processes:
        process.start()
