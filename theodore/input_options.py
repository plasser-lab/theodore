"""
Utilities for reading and writing options from/to an input file.
"""

from __future__ import print_function, division
import sys
from . import error_handler, lib_file


def user_input_py2(inpstr):
    return raw_input(inpstr)


def user_input_py3(inpstr):
    return input(inpstr)


if sys.version_info[0] == 2:
    user_input = user_input_py2
else:
    user_input = user_input_py3


class options:
    """
    Base class for handling input options.
    """
    def __init__(self, ifile):
        self.opt_dict = {}
        self.descr_dict = {} # Collect information for documentation
        self.doc_list = []
        self.ifile = ifile

    def __getitem__(self, option):
        return self.get(option, strict=True)

    def get(self, option, strict=True):
        """
        Return the value of an option.
        """
        self.chk_option(option)

        if strict and self.opt_dict[option] == None:
            raise error_handler.MsgError('Option "%s" not defined in file %s!'%(option, self.ifile))
        else:
            return self.opt_dict[option]

    def __setitem__(self, key, val):
        self.opt_dict[key] = val

    def set_kd(self, key, val, descr=''):
        """
        Set key, value along with a description.
        """
        self.opt_dict[key]   = val
        self.descr_dict[key] = descr
        self.doc_list.append(key)

    def __contains__(self, option):
        """
        Check if an option has been set.

        Raise an error if the option does not even exist.
        """
        self.chk_option(option)

        return self.opt_dict[option] != None

    def has_key(self, option):
        return self.__contains__(option)

    def chk_option(self, option):
        if not option in self.opt_dict:
            raise error_handler.MsgError("Option %s not known!"%option)

    def check_at_lists(self, at_lists, prt_lvl=0):
        """
        Check if an at_lists definition of molecular fragments is useful.
        """
        num_lists = len(at_lists)

        lens = []
        sum_list = []
        for at_list in at_lists:
            sum_list+=at_list
            lens.append(len(at_list))

        numen = len(sum_list)
        maxen = max(sum_list)

        if prt_lvl >= 1:
            print('\nChecking whether the at_lists definition is valid ...')
            if prt_lvl >= 2:
                print('at_lists=', at_lists)
            print('  %i lists with individual numbers of entries:'%(num_lists))
            print(lens)

            print('  %i total entries, with maximal value %i'%(numen,maxen))

        for i in range(1,maxen+1):
            ci = sum_list.count(i)
            if ci!=1:
                print(' WARNING: value %i present %i times in at_lists!'%(i,ci))

    def copy(self, coptions):
        """
        Copy information from a different options instance.
        """
        for key, val in coptions.opt_dict.items():
            self[key] = val

    def doc_info(self):
        """
        Print documentation of keywords.
        Output is in rst format, which can be read directly or via sphinx.
        """
        wt = lib_file.rsttable(3, colw=[20,15,70])
        wt.add_row(['Keyword', 'Default', 'Description'])
        for key in self.doc_list:
            val = self.get(key, strict=False)
            descr = self.descr_dict[key]
            wt.add_row([key, val, descr])
        return wt.ret_table()

class read_options(options):
    """
    General class for handling input options read from file.
    """
    def __init__(self, ifile, check_init=True):
        options.__init__(self, ifile)

        self.set_defaults()
        self.init = self.read_ifile()

        if check_init: self.check_init()

        self.post_process()

    def check_init(self):
        """
        Check if the instance was properly initialized (the file was read).
        """
        if self.init > 0:
            print("\n ERROR: Input file %s not found!"%self.ifile)
            print("  Please create this file using theoinp")
            exit(0)

    def set_defaults(self):
        """
        Set defaults for the options.
        All possible options should appear here.
        -> inherit for specific implementations
        """
        pass

    def read_ifile(self):
        """
        Read the input file self.ifile.
        Key and value are separated by '='.
        Leading and trailing whitespace is removed.
        """
        if self.ifile is None:
            return 0
        try:
            fileh = open(self.ifile, 'r')
        except:
            return 1

        for line in fileh:
            # take out possible comments
            if '#' in line: continue

            words = line.strip().split('=')
            if len(line.strip()) == 0: continue

            if len(words) != 2:
                print(" ERROR: in file %s\n   line cannot be parsed:"%self.ifile)
                print(len(line))
                print(line)
                exit(6)

            key = words[0].strip()

            if words[1] == '':
                raise error_handler.MsgError('Please specify a value for "%s=" in %s!'%(key, self.ifile))

            val = eval(words[1])

            # every possible option has to be initiliazed in set_defaults to avoid confusion
            if not key in self.opt_dict:
                raise error_handler.MsgError('Unknown option in %s: %s'%(self.ifile, key))

            self.opt_dict[key] = val

        return 0

    def get_def(self, option, default):
        self.chk_option(option)

        if self.opt_dict[option] == None:
            return default
        else:
            return self.opt_dict[option]

    def post_process(self):
        pass

class write_options(options):
    """
    General class for writing options to an input file.
    """
    def __init__(self, ifile):
        options.__init__(self, ifile)

        self.ostr = ''

    def read_str(self, title, key, *args, **kwargs):
        """
        Read a string from input.
        """
        titlek = "%s (%s):"%(title, key)

        val = self.ret_str(titlek, *args, **kwargs)

        self.write_option(key, val)

    def ret_str(self, title, default='', autocomp=False):
        # readline, which is used for auto completion,
        #   creates weird output of the form [?1034h
        #   it should only be imported here
        import readline

        print()
        print(title)

        acstr = ' (autocomplete enabled)' if autocomp else ''
        inpstr = 'Choice%s: '%acstr
        if not default=='': inpstr += '[%s] '%default

        if autocomp:
            readline.set_completer_delims(' \t\n;')
            readline.parse_and_bind("tab: complete")    # activate autocomplete
        val = user_input(inpstr)
        readline.parse_and_bind("tab: ")            # deactivate autocomplete

        if val=='': val = default

        return val

    def read_float(self, title, key, default=1.111):
        """
        Read a float from input.
        """
        titlek = "%s (%s):"%(title, key)

        val = self.ret_float(titlek, default)

        self.write_option(key, val)

    def ret_float(self, title, default=1.111):
        print()
        print(title)

        inpstr = 'Choice: '
        if not default==1.111: inpstr += '[%f] '%default

        sval = user_input(inpstr)
        if sval=='':
            val = default
        else:
            val = float(sval)

        return val

    def read_int(self, title, key, idef=-1):
        """
        Read a string from input.
        """
        titlek = "%s (%s):"%(title, key)

        val = self.ret_int(titlek, idef)

        self.write_option(key, val)

    def ret_int(self, title, idef=-1):
        print()
        print(title)

        return self.inp_int(idef)

    def inp_int(self, idef=-1):
        inpstr = 'Choice: '
        if not idef==-1:
            inpstr += '[%i] '%idef

        retval = idef
        while True:
            try:
                retval = int(user_input(inpstr))
            except:
                if retval==-1:
                    print("Please enter an integer number!")
            if retval!=-1: break

        return retval

    def read_yn(self, title, key, default=False):
        """
        Read Boolean from input.
        """
        titlek = "%s (%s):"%(title, key)

        val = self.ret_yn(titlek, default)

        self.write_option(key, val)

        return val

    def ret_yn(self, question, default=False):
        """
        Ask a yes/no question and return True or False.
        """
        print()
        print(question)

        inpstr = 'Choice (y/n): '
        if default:
            inpstr += '[y] '
        else:
            inpstr += '[n] '

        answer = user_input(inpstr)

        if default:
            return not 'n' in answer.lower()
        else:
            return 'y' in answer.lower()

    def choose_list(self, title, key, opt_expl, default=''):
        """
        Choose an option from a list containing options and explanations.
        """
        titlek = "%s (%s):"%(title, key)

        expl = ["%10s - %s"%(opt, expl) for opt, expl in opt_expl]

        idef = -1
        for ioe, oe in enumerate(opt_expl):
            if oe[0] == default:
                idef = ioe + 1
        ichoice = self.ret_choose_list(titlek, expl, idef)

        val = opt_expl[ichoice-1][0]

        self.write_option(key, val)

    def ret_choose_list(self, title, expl, idef=-1):
        """
        Choose an option from a list containing explanations and return the answer.
        """
        print()
        print(title)

        self.print_list(expl)

        return self.inp_int(idef)

    def print_list(self, plist):
        """
        Print an indexed list to screen.
        """
        iopt = 0
        for p in plist:
            iopt += 1
            print("  [%2i] %s"%(iopt, p))

    def write_list(self, key, wlist, lformat="%i"):
        # write_option can be called directly
        self.write_option(key, wlist)

    def write_option(self, key, val):
        self[key] = val

        if type(val) is str:
            self.ostr += "%s='%s'\n"%(key, str(val))
        else:
            self.ostr += "%s=%s\n"%(key, str(val))

    def flush(self, lvprt=0, choose_file=False):
        if choose_file:
            act_ifile = self.ret_str('Name of input file', self.ifile)
        else:
            act_ifile = self.ifile

        fileh = open(act_ifile, 'w')
        fileh.write(self.ostr)
        fileh.close()
        if lvprt==1:
            print('Finished: File %s written.'%act_ifile)

class dens_ana_options(read_options):
    """
    Input options for density analysis.
    """
    def set_defaults(self):
        # General options
        self.set_kd('lvprt', 1, 'Print level')

        # Read options
        self.set_kd('mo_file', None, 'MO-coefficient file (Molden format)')
        self.set_kd('rtype', None, 'Third party program and calculation type')
        self.set_kd('rfile', None, 'Main output file of the third party program')
        self.set_kd('rfile2', None, 'Second output file (if applicable)')
        self.set_kd('ana_files', [], 'List of files to analyze')
        self.set_kd('read_binary', False, 'Read information from a binary file (if applicable)')
        self.set_kd('read_libwfa', False, 'Switch to libwfa output (applicable for qctddft, rassi)')
        self['s_or_t'] = None, # State or transition density matrix analysis (internal only)
        self.set_kd('ignore_irreps', [],  'Ignore irreps in the MO file')
        self.set_kd('min_bf', (), 'Min. contrib. of a basis function type in the MO file, e.g. (2, 0.5)')
        self.set_kd('rd_ene', False, 'Interpret energies as occupations in the NO files')
        self.set_kd('occ_fac', 1., 'Multiply NO occpuations by this factor')
        self.set_kd('unrestricted', False,  'Read unrestricted orbitals')
        self.set_kd('ana_states', [], 'Analyze only a set of states (list starts with 1)')

        # atoms
        self.set_kd('at_lists', None, 'Fragment definition for CT number analysis')
        self.set_kd('frag_lists', None, 'List of fragments')
        self.set_kd('coor_file', None, 'File with coordinates')
        self.set_kd('coor_format', None, 'Format of coordinate file as defined by Open Babel')

        # Output options
        self.set_kd('output_file', "ana_summ.txt", 'Main output file')
        self.set_kd('jmol_orbitals', True, 'Export orbitals (NTOs or NDOs) as Jmol script')
        self.set_kd('molden_orbitals', False, 'Create a Molden file for each set of NTOs or NDOs')
        self.set_kd('min_occ', 0.05, 'Minimal occupation for orbital print out')
        self.set_kd('alphabeta', False, 'Use alpha/beta rather than neg./pos. to code for hole/electron?')
        self.set_kd('mcfmt', '% 10E', 'Output format for Molden coefficients')
        self.set_kd('output_prec', (7,3), 'Number of total and decimal digits for print-out of the summary')
        self.set_kd('print_sorted', True, 'Print the final output sorted by energies')

        # tden analysis
        self.set_kd('Om_formula', 1, 'How to compute Omega: 0 - Mulliken (simple), 1 - Mulliken, 2 - Lowdin')
        self.set_kd('prop_list', [], 'List of properties for final print out.')
        self.set_kd('print_OmAt', False, 'Print the OmAt to .npy file and use for automatic restart')
        self.set_kd('print_OmFrag', True, 'Print out the fragment Omega matrix to an ASCII file')
        self.set_kd('eh_pop', 1, 'Print e/h populations: 1 - for fragments, 2 - also for atoms, 3 - bfs')
        self.set_kd('comp_ntos', True, 'Compute natural transition orbitals')
        self.set_kd('comp_dntos', False, 'Compute the domain NTOs')
        self.set_kd('dnto_frags', [], 'Compute DNTOs only for these fragments')
        self.set_kd('ref_state', (1,1), 'Which state to use as reference (irrep, state)')

        # sden analysis
        self.set_kd('pop_ana', True, 'Perform a population analysis')
        self.set_kd('unpaired_ana', True, 'Perform analysis of unpaired electrons')
        self.set_kd('NO_ana', True, 'Compute natural orbitals')
        self.set_kd('AD_ana', True, 'Perform attachment/detachment analysis and compute the NDOs')
        self.set_kd('BO_ana', True, 'Bond order analysis')
        self.set_kd('min_BO', 0.5,  'Minimal bond order to print')
        self.set_kd('mo_pop_type', -1, 'Pop. ana. of MOs: 1 - for atoms, 2 - for basis function types')

        # options for orbkit
        self.set_kd('cube_orbitals', False, 'Output orbitals as cube files (requires orbkit)?')
        self.set_kd('vmd_ntos', False, 'VMD network for NTOs')
        self.set_kd('vmd_ntos_iv', 0.01, 'isovalue')
        self.set_kd('comp_p_h_dens', False, 'Electron/hole densities as cube files')
        self.set_kd('vmd_ph_dens', False)
        self.set_kd('vmd_ph_dens_iv', 0.01)
        self.set_kd('comp_rho0n', False, 'Transition densities as cube files')
        self.set_kd('vmd_rho0n', False)
        self.set_kd('vmd_rho0n_iv', 0.01)
        self.set_kd('comp_rho', False, 'Densities and unpaired densities as cube files')
        self.set_kd('numproc', 1, 'Num. processors for orbkit')
        self.set_kd('orbkit_extend', 4.0, 'Extension of cube files around molecule')
        self.set_kd('orbkit_step',   0.4, 'Step size in cube file')
        self.set_kd('comp_dnto_dens', 0, 'Cube files for DNTO densities 0 - none, 1 - hole, 2 - elec., 3 - both')
        self.set_kd('fchk_dnto_dens', 0, 'Print DNTO densities to the fchk file (0-3)')
        self.set_kd('normalize_dnto_dens', False, 'Normalize the DNTO densities for each fragment')

        # Additional information
        # irrep labels for output
        self['irrep_labels'] = ['D1', 'D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8']
        self.set_kd('ncore', {}, 'Dictionary: number of frozen core orbitals per irrep')

        # Internal switches (these should not be set in the input file)
        self['spin'] = 0 # 0 - RHF orbitals, 1 - alpha spin, -1 - beta spin

        # Program specific options
        self.set_kd('TDA', False, 'Use of TDA rather than full TDDFT (Q-Chem)')


class tden_ana_options(dens_ana_options):
    """
    Input options for transition density analysis.
    """
    def set_defaults(self):
        dens_ana_options.set_defaults(self)
        # Read options
        self['s_or_t'] = 't'

        # Output options
        self['output_file']   = "tden_summ.txt"
        self['prop_list'] = ['Om', 'POS', 'PR', 'CT', 'COH', 'CTnt']

class sden_ana_options(dens_ana_options):
    """
    Input options for state density analysis.
    """
    def set_defaults(self):
        dens_ana_options.set_defaults(self)
        # Read options
        self['s_or_t'] = 's'

        # Output options
        self['output_file']   = "sden_summ.txt"
        self['prop_list'] = ['nu', 'nunl', 'y0', 'y1', 'eta', 'p']

class libwfa_parse_options(dens_ana_options):
    """
    Input for parsing libwfa output.
    """
    def set_defaults(self):
        dens_ana_options.set_defaults(self)

        self['output_file']   = "libwfa_summ.txt"
        self['prop_list'] = ['Om_', 'nu', 'nunl', 'PRNTO', 'Z_HE', 'dexc', 'dH-E', 'sigH', 'sigE', 'COV', 'Corr', 'p', 'sigD', 'sigA']
        self['read_libwfa'] = True

class fcd_ana_options(dens_ana_options):
    """
    Input options for fragment charge difference analysis.
    """
    def set_defaults(self):
        dens_ana_options.set_defaults(self)

        self['state_pair'] = None
