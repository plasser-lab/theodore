"""
General tools for using pytest.
"""

import os, shutil, subprocess, difflib
import io
import sys
from contextlib import contextmanager
from .actions import ActionFactory


class StringIO_(io.StringIO):

    def read(self):
        pos = super().tell()
        super().seek(0)
        result = super().read()
        super().seek(pos)
        return result


@contextmanager
def mock_stdout():
    old = sys.stdout
    sys.stdout = StringIO_()
    yield sys.stdout
    sys.stdout = old


@contextmanager
def mock_stdin(string):
    old = sys.stdin
    sys.stdin = StringIO_(string)
    yield sys.stdin
    sys.stdin = old


@contextmanager
def commandline(string):
    old = sys.argv
    sys.argv = string.split()
    yield sys.argv
    sys.argv = old


class pytest_job:
    """
    Run and check job in EXAMPLES directory using pytest.
    """
    def __init__(self, cfile, thresh=1.e-6):
        self.wstring = ''
        self.epath = os.path.dirname(os.path.abspath(cfile))
        self.prep()
        self.thresh = thresh

    def run_standard(self):
        """
        Run tests in standard format.
        """
        os.chdir(self.epath + '/RUN')
        for ifile in sorted(os.listdir('../IN_FILES')):
            shutil.copy("../IN_FILES/"+ifile, ifile)
            col = ifile.split('.')
            dtype, atype = col[0], col[2]
            with commandline(f'theodore analyze_{dtype} -f {ifile}'), mock_stdout() as out:
                ActionFactory.from_commandline()
                outlines = out.read()
                with open(f'analyze_{atype}.out', 'w') as fh:
                    fh.write(outlines)
        self.check()

    def run_util(self, args, stdin='', outf=None):
        """
        This is for running an arbitrary utility script.
        Generate the input for stdin as, e.g.:
            tee plot_omfrag.in | theodore plot_omfrag
        After this, self.check() has to be run.
        """
        os.chdir(f'{self.epath}/RUN')
        with commandline(f'theodore {args}'), mock_stdin(stdin), mock_stdout() as out:
            ActionFactory.from_commandline()
            if not outf is None:
                open(outf, 'w').write(out.read())

    def prep(self):
        """Create `RUN` dir"""
        os.chdir(self.epath)
        if os.path.exists('RUN'):
            shutil.rmtree('RUN')
        shutil.copytree('QC_FILES', 'RUN')
        os.chdir(self.epath + '/RUN')

    def check(self):
        """ Check if there are any differences.  """
        os.chdir(self.epath + '/RUN')
        print("Checking primary output files")
        for rfile in os.listdir('../REF_FILES'):
            print("  -> " + rfile)
            self.file_diff('../REF_FILES/'+rfile, rfile)

        if len(self.wstring) > 0:
            raise pytestDiffError(self.epath + '\n\n' + self.wstring)

    def file_diff(self, reff, runf):
        """
        Check if the files are different.
         line-by-line treatment for txt files, otherwise general diff.
        """
        wstring = "\n"

        ref = open(reff).readlines()
        run = open(runf).readlines()
#        suffix = runf.split('.')[-1]

        # These files should match line-by-line
        if '.txt' in runf:
            for iline, line in enumerate(ref):
                if not line == run[iline]:
                    print(' *** Running numerical diff ***')
                    if self.num_diff(line, run[iline]) > 0:
                        self.wstring += 'Numerical threshold exceeded \n\n'
                        diffl = list(difflib.unified_diff(ref, run, fromfile=reff, tofile=runf))
                        wstring = "\n"
                        for line in diffl:
                            self.wstring += line
                        return
        # Use a general diff here
        else:
            diffl = list(difflib.unified_diff(self.diff_ignore(ref), self.diff_ignore(run), fromfile=reff, tofile=runf))
            if len(diffl) > 0:
                wstring = "\n"
                for line in diffl:
                    self.wstring += line

    def diff_ignore(self, dlist):
        outl = []
        iglist = ["TheoDORE", "time", "rbkit", "openbabel", "capabilities", "cclib", "Menger"]
        for line in dlist:
            for ig in iglist:
                if ig in line:
                    break
                elif line.strip() == "":
                    break
            else:
                outl.append(line)
        return outl

    def num_diff(self, rline, line):
        """
        Run numerical diff for line.
        """
        ierr = 0

        rwords = rline.split()
        words  = line.split()

        if not len(rwords) == len(words):
            return 2

        for i, rword in enumerate(rwords):
            try:
                rval = float(rword)
            except ValueError:
                continue
            val = float(words[i])
            if abs(rval-val) > self.thresh:
                ierr += 1

        return ierr

class pytestDiffError(Exception):
    def __init__(self, errmsg):
        self.errmsg = errmsg

    def __str__(self):
        return "\n\n *** pytest detected difference to reference results: ***\n\n %s"%self.errmsg
