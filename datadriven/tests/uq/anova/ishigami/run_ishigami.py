'''
Created on Sep 19, 2016

@author: franzefn
'''
from argparse import ArgumentParser
from multiprocessing.process import Process
from uq.anova.ishigami.test_ishigami import run_ishigami_pce, \
    run_ishigami_sg

if __name__ == '__main__':
    parallel = False

    parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
    parser.add_argument('--parallel', default=False, action='store_true', help='run in parallel')
    args = parser.parse_args()
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
    for gridType in scenarions_sg['gridType']:
        for level in scenarions_sg['level']:
            for refinement in scenarions_sg['refinement']:
                if args.parallel:
                    myargs = (gridType, level, 240, False, refinement, True)
                    processes.append(Process(target=run_ishigami_sg, args=myargs))
                else:
                    run_ishigami_sg(gridType,
                                    level,
                                    200,
                                    False,
                                    refinement,
                                    True)

    for sampler in scenarions_pce['sampler']:
        for degree_1d in scenarions_pce['degree']:
            if args.parallel:
                myargs = (sampler, degree_1d, True)
                processes.append(Process(target=run_ishigami_pce, args=myargs))
            else:
                run_ishigami_pce(sampler, degree_1d, True)

#     # run applications in parallel if there are any available
#     for process in processes:
#         process.start()
