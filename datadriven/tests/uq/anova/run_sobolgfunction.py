'''
Created on Sep 19, 2016

@author: franzefn
'''
from argparse import ArgumentParser
from multiprocessing.process import Process
from uq.anova.test_sobolgfunction import run_sobol_g_function_pce, \
    run_sobol_g_function_sg

if __name__ == '__main__':
    parallel = False

    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--parallel', default=False, action='store_true', help='run in parallel')
    args = parser.parse_args()
    scenarions_sg = {'full': {'gridType': ["linear",
                                           "modlinear",
                                           "modPolyClenshawCurtis",
                                           "modPoly"],
                              'level': [2, 3],
                              'refinement': ['var', 'simple',
                                             'exp', 'squared']},
                     'reduced': {'gridType': ["linear",
                                              "modlinear",
                                              "modPolyClenshawCurtis",
                                              "modPoly"],
                                 'level': [2],
                                 'refinement': ['var', 'simple',
                                                'exp', 'squared']}}
        
    scenarions_pce= {'full': {'sampler': ['fekete', 'leja'],
                              'degree': [2, 3, 4]},
                     'reduced': {'sampler': ['fekete', 'leja'],
                                 'degree': [2, 3, 4, 5, 6, 7, 8]}}

    processes = []
    for model, surrogate in scenarions_sg.items():
        for gridType in surrogate['gridType']:
            for level in surrogate['level']:
                for refinement in surrogate['refinement']:
                    if args.parallel:
                        myargs = (model == "full",
                                  gridType, level, 200, False, refinement)
                        processes.append(Process(target=run_sobol_g_function_sg, args=myargs))
                    else:
                        run_sobol_g_function_sg(model == "full",
                                                gridType,
                                                level,
                                                200,
                                                False,
                                                refinement)

    for model, surrogate in scenarions_pce.items():
        for sampler in surrogate['sampler']:
            for degree_1d in surrogate['degree']:
                if args.parallel:
                    myargs = (model == "full", sampler, degree_1d)
                    processes.append(Process(target=run_sobol_g_function_pce, args=myargs))
                else:
                    run_sobol_g_function_pce(model == "full", sampler, degree_1d)

#     # run applications in parallel if there are any available
#     for process in processes:
#         process.start()
