"""
Distribution for shared-disk systems. Code/Working directory should live on a
Filesystem visible from any host (e.g. NFS-home)

The UQSetting need to have an attribute 'setupCommand' which is a string
with python code. The code should load a UQSetting with proper
Simulation set.
"""
import os
from pysgpp.extensions.datadriven.uq.sampler import Sample
import json

from time import sleep
import random

# username
username = "franzefn"
# username = "leichtdk"

# hostname: number of workers
hosts = {
#          'localhost': {'cores': 4,
#                        'free_cores': 2,
#                        'scratch': '/data2/scratch/Fabian/localhost'},
#         'helium': {'cores': 36,  # actually 42 cores, but I don't use all of them
#                    'free_cores' : 4,
#                    'scratch': '/data2/scratch/%s/helium' % username},
         'kepler': {'cores': 24,
                    'free_cores' : 24,
                    'scratch': '/scratch/sgs/%s' % username}
         #'vgpu1':  {'cores': 12,
         #           'free_cores' : 12,
         #           'scratch': '/home/%s/data2/scratch/%s/vgpu1' % (username, username)}
         # 'vgpu2':  12,
         # 'pcsgs05': 3,
         # 'pcsgs02': 4,
         # 'pcsgs03': 2,
         # 'pcsgs04': 4,
         }


def choose_host():
    """
    selects a host randomly. Every host is replicated
    with respect to the number of cores it has.
    """
    
    """
    allslots = []
    for h, props in hosts.iteritems():
        #allslots = allslots + [h] * props['cores']
        allslots = allslots + [h] * props['free_cores']
    
    while True:
        if not allslots == []:
            host = random.choice(allslots)
            # hosts[host] = hosts[host] - 1
            hosts[host]['free_cores'] = hosts[host]['free_cores'] - 1
            sleep(2)  # do not kill the server
            return host
        else:
            sleep(2)  # wait for cores, that may become free
    """
    allslots = []
    sleep_time = 0.5  # seconds
    for h, props in hosts.iteritems():
        allslots = allslots + [h] * props['free_cores']
    if not allslots == []:
        host = random.choice(allslots)
        hosts[host]['free_cores'] = hosts[host]['free_cores'] - 1
        sleep(sleep_time)  # do not kill the server
        return host
    else:
        sleep(sleep_time)  # wait for cores, that may become free
        return ""


def free_host(host):
    # probably do some load balancing here.
    #pass
    # hosts[host] = hosts[host] + 1
    hosts[host]['free_cores'] = hosts[host]['free_cores'] + 1


def do_dist_run(host, uqsetting, outfile, samplelist, starti, scratch_path, sample_number):
    """
    Calculate the given samples on a remote machine.
    exec's, call after fork'ing. (Does not return)
    @param host: string, name of the host
    @param uqsetting: UQSetting
    @param outfile: string, location where the UQSetting of the result is
                     stored
    @param samplelist: Samples, which samples should be run
    @param starti: int
    """
    # compute string represenation of samples using json
    samples = repr([s.toJson() for s in samplelist])
    # mask the quotation marks with escape sequence
    samples = samples.replace('"', '\\' + '"')
    samples = samples.replace("'", '"')

    infile = repr(uqsetting.getFilename())
    out = repr(outfile)
    
    hostname = repr(host)
    
    scratch = repr(scratch_path)
    
    sam_no = repr(sample_number)

    # the username is different, so the path has to be modified
    pwd = os.getcwd().split('/')
    pwd[2] = username
    pwd = repr("/".join(pwd))

    setup = repr(uqsetting.setupCommand)

    # replace single quotation marks with double ones ->
    # enables right interpretation of strings when forwarding
    # the command to the remote worker
    pwd = pwd.replace("'", '"')
    setup = setup.replace("'", '"')
    infile = infile.replace("'", '"')
    out = out.replace("'", '"')
    hostname = hostname.replace("'", '"')
    scratch = scratch.replace("'", '"')
    sam_no = sam_no.replace("'", '"')

    # python command
    data = "from pysgpp.extensions.datadriven.uq.uq_setting.remote_worker import *; \
    pwd = %s; \
    setup = %s; \
    infile = %s; \
    out = %s; \
    samples = %s; \
    starti = %s; \
    hostname = %s; \
    scratch = %s; \
    sample_number = %s; \
    dist_main(pwd, setup, infile, out, samples, starti, hostname, scratch, sample_number)" \
    % (pwd, setup, infile, out, samples, starti, hostname, scratch, sam_no)

    # now connect per ssh and run the code on host
    # -> the method dist_main is called
    print "=====> new run on %s" % host
    print "It will be executed: %s" % data

    # debugging
    #pdb.set_trace()

    os.execlp("ssh", "ssh", host, "python2 -c '%s'" % data)

    # this will never return


def dist_main(pwd, setup, uq_in, uq_out, samples, starti, host, scratch_path, sample_number):
    """
    Main program for the worker.
    @param pwd: string, working directory
    @param setup: string, python command that builds new UQSetting
    @param uq_in:
    @param uq_out:
    @param samples:
    @param starti:
    """
    os.chdir(pwd)
    print "Worker ", sample_number, " in %s" % pwd

    # Parses the list of sample representations to a list of "real" samples
    # print "=====> samplelist with strings: ", samples[0]  # testing
    samplelist = [json.loads(s) for s in samples]
    # print "=====> samplelist with jsonObjects: ", samplelist[0]  # testing
    samplelist = [Sample.fromJson(jsonObject) for jsonObject in samplelist]

    # build environment for setup
    var = {}
    var['filename'] = uq_in
    var['simResultsToPath_in'] = scratch_path

    # execute setup command -> build uqSetting and saves it
    # to the environment specified by var
    exec setup in var
    uq = var['uq']

    uq.setFilename(uq_out)

    uq.setSimulationStats({})
    uq.setPreprocessorStats({})
    uq.setPreprocessorStatsReverse({})
    uq.setPostprocessorStats({})

    # run the samples
    #uq.runSamples_withoutDistribution(samplelist, starti)
    uq.runSamples_withoutDistribution(samplelist)
    uq.writeToFile()

    # testing
    print "=====> Worker %d finished." % os.getpid()
