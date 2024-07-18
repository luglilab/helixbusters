import sys, traceback
import optparse
from core import *


parser = optparse.OptionParser(usage='', version='1.0')

parser.add_option('-s', action="store", dest="samplesheet", default="",
                  help='ciao ciao')


options, args = parser.parse_args()

if __name__ == '__main__':
    DictInfo = dict()

    run = Helixbusters(samplesheet = options.samplesheet)
