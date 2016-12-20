import subprocess
from argparse import ArgumentParser

parser = ArgumentParser(description='Get a program and run it with input', version='%(prog)s 1.0')
parser.add_argument('--surrogate', default="sg", type=str, help="define which surrogate model should be used (sg, pce)")

scriptname = "run_atan.py --surrogate=%s" % args.surrogate
proc = subprocess.Popen(["pkill", "-f", scriptname], stdout=subprocess.PIPE)
proc.wait()
