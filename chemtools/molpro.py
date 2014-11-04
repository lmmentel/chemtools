from code import Code
from subprocess import Popen, PIPE
import os
import re
import sys

class Molpro(Code):
    '''
    Generic class holding a Molpro object.
    '''

    def __init__(self, **kwargs):
        super(Molpro, self).__init__(**kwargs)

    def write_input(self, inpfile=None, core=None, bs=None, inpdata=None, mol=None):

        if isinstance(bs, list):
            basstr = "".join(x.write_molpro() for x in bs)
        else:
            basstr = bs.write_molpro()
        inpdata = re.sub('geometry', mol.molpro_rep(), inpdata, flags=re.I)
        inpdata = re.sub('basis', "basis={\n"+basstr+"\n}\n", inpdata, flags=re.I)
        if core:
            inpdata = re.sub("core","core,{0:s}\n".format(",".join([str(x) for x in core])), inpdata, flags=re.I)
        else:
            inpdata = re.sub("core","", inpdata, flags=re.I)

        with open(inpfile, 'w') as inp:
            inp.write(inpdata)

    def run(self, inpfile):
        '''
        Run a single molpro job interactively - without submitting to the queue.
        '''

        if "-o" in self.runopts:
            outfile = self.runopts[self.runopts.index("-o") + 1]
        else:
            outfile = os.path.splitext(inpfile)[0] + ".out"
        errfile = os.path.splitext(outfile)[0] + ".err"
        opts = []
        opts.extend([self.execpath, inpfile] + self.runopts)

        process = Popen(opts, stdout=PIPE, stderr=PIPE)
        out, err = process.communicate()
        ferr = open(errfile, 'w')
        ferr.write(out)
        ferr.write("{0:s}\n{1:^80s}\n{0:s}\n".format("="*80, "Error messages:"))
        ferr.write(err)
        ferr.close()

        return outfile

    def run_multiple(self, inputs):
        '''
        Run a single molpro job interactively - without submitting to the queue.
        '''

        procs = []
        outputs = [os.path.splitext(inp)[0] + ".out" for inp in inputs]
        for inpfile, outfile in zip(inputs, outputs):
            opts = []
            opts.extend([self.execpath, inpfile] + self.runopts)
            out = open(outfile, 'w')
            process = Popen(opts, stdout=out, stderr=out)
            out.close()
            procs.append(process)

        for p in procs: p.wait()

        return outputs

    def parse(self, output, method, objective, regexp=None):
        '''
        Parser molpro output file to get the objective.
        '''

        parser = MolproOutputParser(output)

        if objective == "total energy":
            if method == "hf":
                return parser.get_hf_total_energy()
            elif method == "cisd":
                return parser.get_cisd_total_energy()
        elif objective == "correlation energy":
                return parser.get_cisd_total_energy() - parser.get_hf_total_energy()
        elif objective == "core energy":
            if method == "cisd":
                return parser.get_cisd_total_energy()
        elif objective == "regexp":
            return parser.get_variable(regexp)
        else:
            sys.exit("<parse>: unknown objective {0:s}".format(objective))

    def accomplished(self, outfile=None):
        '''
        Return True if Molpro job finished without errors.
        '''

        if outfile is not None:
            parser = MolproOutputParser(outfile)
        else:
            parser = MolproOutputParser(self.outfile)
        return parser.terminatedOK()

    def __repr__(self):
        return "<Molpro(\n\tname='{n}',\n\texecpath='{e}',\n\trunopts='{r}')>".format(n=self.name, e=self.execpath, r=self.runopts) 

class MolproOutputParser(object):

    '''Class for parsing molro output files'''

    def __init__(self, out=None):

        self.output = out
        self.outexists()

    def outexists(self):

        '''Check if the out file exists.'''

        if os.path.exists(self.output):
            return True
        else:
            sys.exit("Molpro out file: {0:s} doesn't exist in {1:s}".format(
                     self.output, os.getcwd()))

    def get_hf_total_energy(self):

        '''Return the total HF energy.'''

        with open(self.output, 'r') as out:
            data = out.read()

        hfre = re.compile(r'!RHF STATE \d+\.\d+ Energy\s+(?P<energy>\-?\d+\.\d+)', flags=re.M)
        match = hfre.search(data)
        if match:
            return float(match.group("energy"))

    def get_mp2_total_energy(self):

        '''Return the total MP2 energy.'''

        with open(self.output, 'r') as out:
            data = out.read()

        mpre = re.compile(r'!MP2 total energy\s+(?P<energy>\-?\d+\.\d+)', flags=re.M)
        match = mpre.search(data)
        if match:
            return float(match.group("energy"))

    def get_ccsd_total_energy(self):

        '''Return the total CCSD energy.'''

        with open(self.output, 'r') as out:
            data = out.read()

        ccre = re.compile(r'!CCSD total energy\s+(?P<energy>\-?\d+\.\d+)', flags=re.M)
        match = ccre.search(data)
        if match:
            return float(match.group("energy"))

    def get_ccsdt_total_energy(self):

        '''Return the total CCSD(T) energy.'''

        with open(self.output, 'r') as out:
            data = out.read()

        ccre = re.compile(r'!CCSD\(T\) total energy\s+(?P<energy>\-?\d+\.\d+)', flags=re.M)
        match = ccre.search(data)
        if match:
            return float(match.group("energy"))

    def get_cisd_total_energy(self):

        '''Return the total CISD energy.'''

        with open(self.output, 'r') as out:
            data = out.read()

        cire = re.compile(r'!(RHF-R)?CISD\s+(total\s+)?energy\s+(?P<energy>\-?\d+\.\d+)', flags=re.M)
        match = cire.search(data)
        if match:
            return float(match.group("energy"))

    def get_fci_total_energy(self):

        '''Return the total HF energy.'''

        with open(self.output, 'r') as out:
            data = out.read()

        fcire = re.compile(r'!FCI STATE \d+\.\d+ Energy\s+(?P<energy>\-?\d+\.\d+)', flags=re.M)
        match = fcire.search(data)
        if match:
            return float(match.group("energy"))

    def get_variable(self, rawstring):
        with open(self.output, 'r') as out:
            data = out.read()

        genre = re.compile(rawstring, flags=re.M)
        match = genre.search(data)
        if match:
            return float(match.group(1))

    def terminatedOK(self):

        '''Check if the job terminated succesfully.'''

        with open(self.output, 'r') as out:
            data = out.read()

        errorre = re.compile(r'\s*error', flags=re.I)

        match = errorre.search(data)
        if match:
            return False
        else:
            return True

def parse_basis(string):
    '''
    Parse basis set from a string in Molpro format.
    '''

    bas_re = re.compile(r'basis\s*=\s*\{(.*?)\}', flags=re.DOTALL|re.IGNORECASE)

    m = bas_re.search(string)
    if m:
        lines = m.group(1).split("\n")
    else:
        raise ValueError("basis string not found")

    start = []
    for i, line in enumerate(lines):
        if line.split(",")[0].lower() in ["s","p", "d", "f", "g", "h", "i"]:
            start.append(i)
    if len(start) == 0:
        return None

    startstop = []
    for i in range(len(start)-1):
        startstop.append((start[i], start[i+1]))
    startstop.append((start[-1], len(lines)))

    bs = {}
    for i in startstop:
        at_symbol, shell = parse_shell(lines[i[0]], lines[i[0]+1:i[1]])
        if at_symbol in bs.keys():
            bs[at_symbol] = dict(list(bs[at_symbol].items()) + list(shell.items()))
        else:
            bs[at_symbol] = shell
    return bs

def parse_shell(expsline, coeffs):
    '''
    Parse functions of one shell in molpro format.
    '''

    fs = {}

    shell  = expsline.split(",")[0]
    at_symbol = expsline.split(",")[1].strip().capitalize()
    exps   = [float(x) for x in expsline.rstrip(";").split(",")[2:]]

    fs[shell.lower()] = {'exponents' : exps, 'contractedfs' : []}
    for line in coeffs:
        lsp = line.rstrip(";").split(",")
        if lsp[0] == "c":
            i, j = [int(x) for x in lsp[1].split(".")]
            coeffs = [float(x) for x in lsp[2:]]
            fs[shell.lower()]['contractedfs'].append({'indices' : list(range(i-1, j)), 'coefficients' : coeffs})
    return at_symbol, fs

def parse_ecp(ecpstring):

    ecp_re = re.compile(r'\!\s*Effective core Potentials.*-{25}\s*\n(.*?)\n\s*\n', flags=re.DOTALL)

    lines = ecpstring.split("\n")

    start = []
    for i, line in enumerate(lines):
        if line.split(",")[0].lower() == 'ecp':
            start.append(i)

    if len(start) == 0:
        return None

    startstop = []
    for i in range(len(start)-1):
        startstop.append((start[i], start[i+1]))
    startstop.append((start[-1], len(lines)))

    ecp = {}
    for i in startstop:
        ecp = dict(list(ecp.items()) + list(parse_coeffs(lines[i[0] : i[1]]).items()))
    return ecp

def parse_coeffs(lines):

    firstl = lines[0].replace(';', '').split(',')
    element = firstl[1].strip().capitalize()
    nele = firstl[2]
    lmax = firstl[3]

    liter = iter(x for x in lines[1:] if x != '')

    ecp = {element : {"nele" : nele, "lmax" : lmax, "shells" : []}}

    while True:
        try:
            temp = next(liter)
        except StopIteration as err:
            break
        nlmax = int(temp.split(";")[0])
        comment = temp.split(";")[1].replace("!", "")
        tt = {'comment' : comment, 'parameters' : []}
        for i in range(nlmax):
            param = next(liter).replace(";", "").split(",")
            tt['parameters'].append({'m' : float(param[0]), 'gamma' : float(param[1]),'c' : float(param[2])})
        ecp[element]['shells'].append(tt)
    return ecp

def write_molpro_basis(basisset):
    '''
    Write basis set in molpro format

    This little function is quite dirty and would benefit from  rewriting!
    '''

    newd = []

    for atom, atomgroup in groupby(basisset, lambda x: x["atomic"]):
        for shell, shellgroup in groupby(atomgroup, lambda x: x["shell"]):
            exps    = []
            indices = []
            coeffs  = []
            for function in shellgroup:
                for zeta in function["exps"]:
                    if zeta not in exps:
                        exps.append(zeta)
                istart = exps.index(function["exps"][0]) + 1
                istop  = exps.index(function["exps"][-1]) + 1
                indices.append((istart, istop))
                coeffs.append(function["coeffs"])
            newd.append({"atomic" : atom,
                            "shell"  : shell,
                            "exps"   : exps,
                            "indices": indices,
                            "coeffs" : coeffs})
    outstring = "basis={\n"
    for f in sorted(newd, key=itemgetter("atomic", "shell")):
        elem = periodic.element(f["atomic"])
        outstring = outstring + "{0}, {1}, {2}\n".format(
                _shells[f["shell"]].lower(),
                elem.symbol,
                ", ".join([str(x) for x in f["exps"]]))
        for i, item in enumerate(f["coeffs"]):
            outstring = outstring + "{0}, {1}, {2}\n".format(
                "c",
                ".".join([str(x) for x in f["indices"][i]]),
                ", ".join([str(x) for x in item]))
    outstring = outstring + "}\n"
    return outstring
