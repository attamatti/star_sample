#!/usr/bin/python
# pick a random sampling of micrographs from a specific range of defocus measurements for testing autopicking


import sys
import collections
import random
import os

#------- function test if string is a number --------------------------#
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
#-----------------------------------------------------------------------

###---------function: read the star file get the header, labels, and data -------------#######
def read_starfile(f):
    alldata = open(f,'r').readlines()
    labelsdic = {}
    data = []
    header = []
    for i in alldata:
        if '#' in i:
            labelsdic[i.split('#')[0]] = int(i.split('#')[1])-1
        if len(i.split()) > 3:
            data.append(i.split())
        if len(i.split()) < 3:
            header.append(i.strip("\n"))
    return(labelsdic,header,data)
#---------------------------------------------------------------------------------------------#

#------ function: write all of the numbers in the fortran format ---------------------------#
def make_pretty_numbers(dataarray):
    prettyarray = []
    for line in dataarray:
        linestr = ""
        for i in line:
            if is_number(i):
                count = len(i.split('.'))
                if count > 1:
                    i = float(i)
                    if len(str(i).split('.')[0]) > 5:
                        linestr= linestr+"{0:.6e} ".format(i)
                    else:
                        linestr= linestr+"{0:12.6f} ".format(i)
                else:
                    linestr= linestr+"{0: 12d} ".format(int(i))
            else:
                linestr= linestr+"{0} ".format(i)
        prettyarray.append(linestr)
    return prettyarray
#---------------------------------------------------------------------------------------------#

#-----------------------------------------------------------------------------------#
errormsg = "USAGE: rln_sample_starfile.py --i <star file>\nOptional arguments:\n--bin <n>\t\t specifiy the size of the defocus bins, default is 0.5 micron"
class Arg(object):
    _registry = []
    def __init__(self, flag, value, req):
        self._registry.append(self)
        self.flag = flag
        self.value = value
        self.req = req

def make_arg(flag, value, req):
    Argument = Arg(flag, value, req)
    if Argument.req == True:
        if Argument.flag not in sys.argv:
            print(errormsg)
            sys.exit("ERROR: required argument '{0}' is missing".format(Argument.flag))
    if Argument.value == True:
        try:
            test = sys.argv[sys.argv.index(Argument.flag)+1]
        except ValueError:
            if Argument.req == True:
                print(errormsg)
                sys.exit("ERROR: required argument '{0}' is missing".format(Argument.flag))
            elif Argument.req == False:
                return False
        except IndexError:
                print(errormsg)
                sys.exit("ERROR: argument '{0}' requires a value".format(Argument.flag))
        else:
            if Argument.value == True:
                Argument.value = sys.argv[sys.argv.index(Argument.flag)+1]
        
    if Argument.value == False:
        if Argument.flag in sys.argv:
            Argument.value = True
        else:
            Argument.value = False
    return Argument.value
#-----------------------------------------------------------------------------------#

#------function return group number --------------#
def returngroupnumber(defo):
    for i in reversed(groupmaxes):
        if float(defo)/10000 >= i:
            return(int(groupmaxes.index(i)))
#------------------------------------------------#

#----- function get random number
def get_rand(range1,range2):
    n = random.randrange(range1,range2)
    if n not in picked:
        return n
    else:
        get_rand(range1,range2)
#----------#  

file = make_arg('--i',True,True)
binsize = make_arg('--bin',True,False)

if os.path.isfile(file) != True:
    sys.exit('ERROR: File {0} not found'.format(file))

(labels,header,data) = read_starfile(file)

defoci = []
for i in data:
    defoci.append(float(i[labels['_rlnDefocusU ']]))


mind =  min(defoci)/10000
maxd =  max(defoci)/10000
if binsize != False:
    groupbin = float(binsize)
else:
    groupbin = 0.5
numgroups = int(((maxd - mind)/groupbin))+1
groupmaxes = [mind]
for i in range (1,numgroups):
    groupmaxes.append(mind+(groupbin*i))
groupids = {}
for i in range(0,numgroups):
    groupids[i] = (round(mind+(groupbin*i),2),round(mind+groupbin*(i+1),2))


groups = {}
for i in data: 
    newdataline = []
    if returngroupnumber(i[labels['_rlnDefocusU ']]) not in groups:
        groups[returngroupnumber(i[labels['_rlnDefocusU ']])] = [i]
    else:
        groups[returngroupnumber(i[labels['_rlnDefocusU ']])].append(i)
keys = range(0,len(groups))
print('   Defocus\tMicrographs # = 10')
zeros = []
for i in keys:
    try:
        count = int(len(groups[i])/10)
    except KeyError:
        count = 'X'
    if count == 0:
        count  = 1
    if count == 'X':
        count = 0
        zeros.append(str(i))
    print('{0}) {1}-{2}\t{3}'.format(i,groupids[int(i)][0],groupids[int(i)][1],count*'#'))

nchosen = raw_input('which groups to use (comma separated): ')
chocheck = nchosen.split(',')
for i in chocheck:
    if int(i) not in range(0,len(keys)+2) and i != ',':
        sys.exit('ERROR: Invalid entry : Out of range')

nmicrographs = int(raw_input('how many from each group: '))

finaldata = []
for i in nchosen.split(','):
    picked = []
    if i in zeros:
        sys.exit('ERROR: No micrographs in group {0}'.format(i))
    if nmicrographs > len(groups[int(i)]):
        sys.exit('ERROR: not enough micrographs in group {0}'.format(i))
    for j in range(0,nmicrographs):
        if len(groups[int(i)]) == 0:
            break
        n = get_rand(0,len(groups[int(i)]))
        finaldata.append(groups[int(i)][n])
        groups[int(i)].remove(groups[int(i)][n])

output = open('random_selection.star','w')

for i in header:
    output.write('{0}\n'.format(i))
final= make_pretty_numbers(finaldata)
for i in final:
    output.write('{0}\n'.format(i))
