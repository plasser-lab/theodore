"""
Input generation for TheoDORE runs.
"""
from __future__ import print_function, division
from .actions import Action
import os
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    input_options = importer.lazy_import_as('..input_options', 'input_options')
    lib_struc = importer.lazy_import_as('..lib_struc', 'lib_struc')
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')
    orbkit_interface = importer.lazy_import_as('..orbkit_interface', 'orbkit_interface')


class write_options_theo(input_options.write_options):
    """
    Input generation with suitable defaults.
    """
    # TODO: this could optionally take information from an existing dens_ana.in file

    def choose_rtype(self):
        rdef = ''
        if os.path.exists('qchem.fchk'):
            rdef = 'fchk'
        elif os.path.exists('qchem.out'):
            rdef = 'qcadc'
        elif os.path.exists('geom'):
            rdef = 'colmrci'
        elif os.path.exists('molcas.rasscf.molden'):
            rdef = 'rassi'
        elif os.path.exists('mrci.cidens') or os.path.exists('mrci.tcidens'):
            rdef = 'dftmrci'
        elif os.path.exists('coord'):
            if os.path.exists('auxbasis'):
                rdef = 'ricc2'
            else:
                rdef = 'escf'
        elif os.path.exists('ams.results/adf.rkf'):
            rdef = 'adf'
        elif os.path.exists('EXC.DAT'):
            rdef = 'tddftb'

        self.choose_list(
            'Type of job',
            'rtype',
        [
            ('qcadc', 'Q-Chem ADC (libwfa output)'),
            ('libwfa', 'General libwfa output'),
            ('qctddft', 'Q-Chem TDDFT'),
            ('fchk', 'Q-Chem fchk file'),
            ('colmcscf', 'Columbus MCSCF'),
            ('colmrci', 'Columbus MR-CI (tden analysis)'),
            ('rassi', 'Molcas RASSI'),
            ('nos', 'Read natural orbitals (Molden format) for sden analysis: Columbus, Molcas, ...'),
            ('ricc2', 'Turbomole ricc2'),
            ('escf', 'Turbomole escf'),
            ('terachem', 'Terachem (TDDFT)'),
            ('cclib', 'Use external cclib library: Gaussian, GAMESS, ...'),
            ('orca', 'ORCA TDDFT (using a Molden file and cclib)'),
            ('orcaCAS', 'ORCA CAS (using .molden, .densities and .densitiesinfo)'),
            ('adf', 'ADF (TDDFT)'),
            ('tddftb', 'DFTB+ - TDDFTB'),
            ('dftmrci', 'DFT/MRCI'),
            ('onetep', 'ONETEP')
        ], rdef)

        # set defaults
        #   None means that the option is not applicable
        #   '' means that no default is given

        self['read_libwfa'] = False
        # optionally this could be read from dens_ana.in
        if self['rtype'] == 'qcadc':
            self['rfile'] = 'qchem.out'
            self['mo_file'] = None
            self['coor_file'] = 'qchem.out'
            self['coor_format'] = 'qcout'
            self['read_libwfa'] = True
        elif self['rtype'] == 'libwfa':
            self['rfile'] = None
            self['mo_file'] = None
            self['coor_file'] = 'qchem.out'
            self['coor_format'] = 'qcout'
            self['read_libwfa'] = True
        elif self['rtype'] == 'qctddft':
            self['rfile'] = 'qchem.out'
            self['mo_file'] = 'qchem.mld'
            self['coor_file'] = 'qchem.out'
            self['coor_format'] = 'qcout'
        elif self['rtype'] == 'fchk':
            self['rfile'] = 'qchem.fchk'
            self['mo_file'] = None
            self['coor_file'] = 'qchem.out'
            self['coor_format'] = 'qcout'
        elif self['rtype'] == 'colmrci' or self['rtype'] == 'colmcscf':
            self['rfile'] = None
            self['mo_file'] = 'MOLDEN/molden_mo_mc.sp'
            self['coor_file'] = 'geom'
            self['coor_format'] = 'col'
        elif self['rtype'] == 'rassi':
            self['rfile'] = 'molcas.log'
            self['mo_file'] = 'molcas.rasscf.molden'
            self['coor_file'] = 'geom.xyz'
            self['coor_format'] = 'xyz'
        elif self['rtype'] == 'nos':
            self['rfile'] = None
            self['mo_file'] = ''
            self['coor_file'] = 'geom'
            self['coor_format'] = 'col'
        elif self['rtype'] == 'ricc2':
            self['rfile'] = 'ricc2.out'
            self['mo_file'] = 'molden.input'
            self['coor_file'] = 'coord'
            self['coor_format'] = 'tmol'
        elif self['rtype'] == 'escf':
            self['rfile'] = 'escf.out'
            self['mo_file'] = 'molden.input'
            self['coor_file'] = 'coord'
            self['coor_format'] = 'tmol'
        elif self['rtype'] == 'terachem':
            self['rfile'] = ''
            self['mo_file'] = 'inpfile.molden'
            self['coor_file'] = ''
            self['coor_format'] = ''
        elif self['rtype'] == 'cclib':
            self['rfile'] = ''
            self['mo_file'] = 'orbs.mld'
            self['coor_file'] = ''
            self['coor_format'] = ''
        elif self['rtype'] == 'orca' or self['rtype'] == 'orcaCAS':
            self['rfile'] = 'orca.out'
            self['mo_file'] = 'orca.molden.input'
            self['coor_file'] = ''
            self['coor_format'] = ''
        elif self['rtype'] == 'adf':
            self['rfile'] = 'ams.results/adf.rkf'
            self['mo_file'] = None
            self['coor_file'] = ''
            self['coor_format'] = ''
        elif self['rtype'] == 'tddftb':
            self['rfile'] = 'EXC.DAT'
            self['mo_file'] = 'eigenvec.out'
            self['coor_file'] = 'geom.xyz'
            self['coor_format'] = 'xyz'
        elif self['rtype'] == 'dftmrci':
            self['rfile'] = 'mrci2.out'
            self['rfile2'] = 'mrci.cidens'
            self['mo_file'] = 'orca.molden.input'
            self['coor_file'] = ''
            self['coor_format'] = ''
        elif self['rtype'] == 'onetep':
            self['rfile'] = ''
            self['mo_file'] = None
            self['coor_file'] = ''
            self['coor_format'] = ''
        else:
            self['rfile'] = ''
            self['mo_file'] = ''
            self['coor_file'] = ''
            self['coor_format'] = ''

    def set_read_options(self):
        if 'rfile' in self:
           self.read_str('Main file to read', 'rfile', self['rfile'], autocomp=True)

        # switch for libwfa
        if self['rtype'] in ['qctddft']:
            if self.read_yn('Did you run "state_analysis=True"?', 'read_libwfa', True):
                self['mo_file'] = None
        elif self['rtype'] in ['rassi']:
            if self.read_yn('Did you use &WFA?', 'read_libwfa', True):
                self['mo_file'] = None
        elif self['rtype'] in ['cclib', 'gamess']:
            print('\nNote: If used in connection with ORBKIT it is preferable to have an externally generated Molden file.')
            if not self.ret_yn('Do you have an externally generated Molden file?', False):
                self['mo_file'] = None

        # switch for TDA
        if self['rtype'] in ['qctddft']:
            self.read_yn('Read TDA rather than full TDDFT results?', 'TDA', False)

        # binary files
        if self['rtype'] == 'ricc2':
            self.read_yn('Read binary CCRE0 files?', 'read_binary', False)
            if not self['read_binary']:
                print()
                print(" *** Warning: without read_binary you have to delete the line")
                print("       implicit core=   x virt=    x")
                print("     from the control file before running tm2molden.")
        elif self['rtype'] == 'orca':
            self.read_yn('Read the binary orca.cis file?', 'read_binary', True)
        elif self['rtype'] == 'orcaCAS':
            self.read_yn('Read the binary orca.densities and orca.densitiesinfo files?', 'read_binary', True)

        # switch for DFTMRCI singlet-singlet / singlet-triplet densities
        if self['rtype'] == 'dftmrci':
            self.read_str('Read mrci.cidens or mrci.tcidens file?', 'rfile2', self['rfile2'], autocomp=True)

        # extra read-out
        if self['rtype'] == 'colmrci':
            if self.ret_yn('Parse extensivity corrected MRCI energies?', True):
                self['prop_list'] += ['+DV3', '+P']

    def choose_mo_file(self):
        if 'mo_file' not in self: return

        mo_str = 'MO file (Molden format)\
                      \n -> This file should ideally contain a square invertible coefficient matrix'
        self.read_str(mo_str, 'mo_file', self['mo_file'], True)

        if self['rtype'] == 'nos':
            nodir = self.ret_str('Directory with the NO files:', '.', True)
            print("%s contains the following files:"%nodir)

            lfiles = sorted(os.listdir(nodir))
            self.print_list(lfiles)
            rstr = self.ret_str('Input indices of required files (separated by spaces)\n Start with ground state.')
            noinds = [int(ino)-1 for ino in rstr.replace(',',' ').split()]

            nolist = [os.path.join(nodir,lfiles[noind]) for noind in noinds]
            self.write_option('ana_files', nolist)

            self.read_yn('Intepret energies as orbital occupations (for Q-Chem)', 'rd_ene', False)

    def only_at_lists(self, afile):
        """
        Shortcut for only creating at_lists
        """
        print("Coordinate file (afile): %s"%afile)
        self['coor_file'] = afile
        tstruc = lib_struc.structure()
        self['coor_format'] = tstruc.guess_file_type(self['coor_file'])
        self.make_at_lists(afile)

        print('\nOpenMolcas input:\n')
        print('&WFA\n H5FI = $Project.rassi.h5\n ATLISTS')
        print(' %i'%len(self['at_lists']))
        for li in self['at_lists']:
            wstr = ' '
            for el in li:
                wstr += '%i '%el
            wstr += '*'
            print(wstr)

    def make_at_lists(self, afile=None):
        print("Fragment definition for CT nubmer analysis")
        aexpl = ['Manual input',\
                 'Automatic generation by fragment (using python-openbabel)', \
                 'Automatic generation for transition metal complexes (using python-openbabel)', \
                 'Mixed manual/automatic generation (using python-openbabel)', \
                 'Automatic generation by element (using python-openbabel)', \
                 'Leave empty and fill out later']
        ichoice = self.ret_choose_list('Mode for specifying molecular fragments (at_lists):', aexpl)

        if ichoice==1:
            self['at_lists'] = self.read_at_lists()
        elif ichoice in [2,3,4,5]:
            self['at_lists'] = self.file_at_lists(ichoice, afile)
        elif ichoice==6:
            self['at_lists'] = [[]]
            self.ostr += 'at_lists=\n'
            return

        self.check_at_lists(self['at_lists'], 2)

        self.ostr += 'at_lists=%s\n'%str(self['at_lists'])

    def read_at_lists(self):
        atl_tmp = []

        for ifrag in range(1, 1000):
            rstr = self.ret_str('Input the indices of the atoms belonging to fragment %i:\n(separated by spaces)'%ifrag)
            if rstr == '': break

            atl_tmp.append([int(iat) for iat in rstr.replace(',', ' ').split()])

        return atl_tmp

    def file_at_lists(self, mode, afile=None):
        print("Automatic generation of at_lists partitioning ...")
        if afile is None:
            self.coor_file()

        struc = lib_struc.structure()
        struc.read_file(file_path=self['coor_file'], file_type=self['coor_format'])

        if mode == 2:
            return struc.ret_partition()
        elif mode == 3:
            rstr = self.ret_str('Input the index of the transition metal atom (or indices of the corresponding fragment)')
            inp_list = [int(iat) for iat in rstr.replace(',', ' ').split()]
            return struc.ret_partition(inp_lists=[inp_list])
        elif mode == 4:
            print("\nEnter the manual fragment defitions first.\nLeave empty to switch to automatic mode.")
            inp_lists = self.read_at_lists()
            return struc.ret_partition(inp_lists=inp_lists)
        elif mode == 5:
            return struc.ret_el_partition()
        else:
            raise error_handler.ElseError('mode', mode)

    def coor_file(self):
        tmp = self['coor_file']
        self.read_str('Coordinate file', 'coor_file', self['coor_file'], True)
        if self['coor_file'] != tmp:
            try:
                tstruc = lib_struc.structure()
                self['coor_format'] = tstruc.guess_file_type(self['coor_file'])
            except error_handler.MsgError:
                pass
        self.read_str('Format of coordinate file', 'coor_format', self['coor_format'])

    def set_Om_desc(self):
        """
        Set the list of Omega descriptors to be computed and set the formula.
        """
        if self['read_libwfa']:
            pass
        elif self['rtype'] == 'onetep':
            self.read_int('Formula for Omega matrix computation\n\
   0 - simple, 1 - Mulliken', 'Om_formula', 1)
        else:
            self.read_int('Formula for Omega matrix computation\n\
   0 - simple, 1 - Mulliken, 2 - Lowdin', 'Om_formula', 2)

        Olist = ['Standard set', 'Transition metal complex', 'None']

        ichoice = self.ret_choose_list('Omega descriptors to be computed:', Olist, 1)

        if ichoice==1:
            self['prop_list'] += ['Om', 'POS', 'PR', 'CT', 'COH', 'CTnt']
        elif ichoice==2:
            self['prop_list'] += ['Om', 'POSi','POSf','PR','CT','MC','LC','MLCT','LMCT','LLCT']

    def set_eh_pop(self):
        """
        Print-out of e/h populations.
        """
        eh_list = ['None', 'For fragments', 'For fragments and individual atoms']

        ichoice = self.ret_choose_list('Print-out of electron/hole populations', eh_list, 1)
        self.write_option('eh_pop', ichoice-1)

    def write_prop_list(self):
        if len(self['prop_list']) == 0: return

        self.write_option('prop_list', self['prop_list'])

    def nto_ana(self):
        if self['read_libwfa']: return

        self.read_yn('Perform natural transition orbital (NTO) analysis?', 'comp_ntos', True)
        self.read_yn('Perform analysis of domain NTOs and conditional densities?', 'comp_dntos', False)
        if self['comp_ntos'] or self['comp_dntos']:
            self['prop_list'] += ['PRNTO', 'Z_HE']
            if self['rtype'] in ['adf', 'tddftb', 'onetep']:
                return
            if self['rtype'] == 'ricc2' and self['read_binary'] == True:
                self.read_yn('NTOs in Molden format', 'molden_orbitals', True)
            else:
                self.read_yn('NTOs as Jmol script?', 'jmol_orbitals', True)
                self.read_yn('NTOs in Molden format', 'molden_orbitals', False)
            if self['molden_orbitals']:
                self.read_yn('Use alpha/beta rather then negative/positive to code for hole/particle orbitals?', 'alphabeta', False)
            if orbkit_interface.orbkit_avail:
                self.read_yn('NTOs in Cube file format (requires orbkit)', 'cube_orbitals', False)
                if self['cube_orbitals']:
                    self.read_yn('Create VMD Network for NTOs', 'vmd_ntos', False)
                    if self['vmd_ntos']:
                        self.read_float('Isosurface value for VMD network', 'vmd_ntos_iv', 0.01)
                self.read_yn('Calculation of Particle/Hole density (requires orbkit)?', 'comp_p_h_dens', False)
                if self['comp_p_h_dens']:
                    self.read_yn('Create VMD Network for p/h densities', 'vmd_ph_dens', False)
                    if self['vmd_ph_dens']:
                        self.read_float('Isosurface value for VMD network', 'vmd_ph_dens_iv', 0.01)
                if self['comp_dntos']:
                    self.read_int('Compute conditional densities as cube files?\n 0 - no, 1 - hole, 2 - electron, 3 - both', 'comp_dnto_dens', 0)
                    if self['rtype'] == 'fchk':
                        self.read_int('Write conditional densities to fchk file\n 0 - no, 1 - hole, 2 - electron, 3 - both', 'fchk_dnto_dens', 0)
        else:
            self.ostr += 'comp_ntos=False\n'

    def exciton_ana(self):
        if self.ret_yn('Compute approximate exciton size?', True):
            self['prop_list'] += ['RMSeh']
            if not self['rtype'] in ['orca', 'cclib', 'adf']:
                if not 'mo_file' in self:
                    print("\nMolecular coordinates for exciton analysis:")
                    self.coor_file()

        if self['read_libwfa']:
            if self.ret_yn('Parse exciton information from libwfa analysis?', False):
                self['prop_list'] += ['dexc', 'dH-E', 'sigH', 'sigE', 'COV', 'Corr']
            if self.ret_yn('Parse 1DDM exciton information from libwfa analysis?', False):
                self['prop_list'] += ['sigD', 'sigA']

    def AD_ana(self):
        self.read_yn('Attachment/detachment analysis', 'AD_ana', True)
        if self['AD_ana']:
            self['prop_list'] += ['p']
            if not self['rtype'] in ['adf', 'tddftb', 'onetep']:
                self.read_yn('NDOs as Jmol script?', 'jmol_orbitals', True)
                self.read_yn('NDOs in Molden format?', 'molden_orbitals', False)
                if self['molden_orbitals']:
                    self.read_yn('Use alpha/beta rather then negative/positive to code for det./att. orbitals?', 'alphabeta', False)
        else:
            self.write_option('jmol_orbitals', 'False')

    def BO_ana(self):
        self.read_yn('Mayer bond order and valence analysis?', 'BO_ana', True)

    def comp_rho0n(self):
        if not orbkit_interface.orbkit_avail:
            return
        if self['rtype'] in ['adf', 'tddftb', 'onetep']:
            return
        if self['read_libwfa'] == True:
            return

        self.read_yn('Calculation of transition densities between ground state and excited states (requires orbkit)', 'comp_rho0n', False)
        if self['comp_rho0n']:
            self.read_yn('Create VMD Network for transition densities', 'vmd_rho0n', False)
            if self['vmd_rho0n']:
                self.read_float('Isosurface value for VMD network', 'vmd_rho0n_iv', 0.01)

    def ddm_parse(self):
        if self['read_libwfa']:
            if self.ret_yn('Parse 1DDM exciton information from libwfa job?', True):
                self['prop_list'] += ['dD-A', 'sigD', 'sigA']

    def get_ncore(self, cidrt_file='cidrtmsin'):
        try:
            lines = open(cidrt_file).readlines()
            ncore = self.parse_ncore_from_cidrtmsin(lines)

            self['ncore_dict'] = ncore
            self.ostr += 'ncore=%s\n' % str(ncore)

            print("Auto-detected frozen core:", ncore)
            return

        except Exception as e:
            print("Auto-detection failed, falling back to manual input:", e)

        if self.ret_yn('Were there frozen core orbitals in the calculation?', True):
            self['ncore_dict'] = {}
            print("Please enter the irrep label and number of orbitals (separated by spaces), e.g. b1u 5")
            for iirrep in range(1, 32):
                rstr = self.ret_str('Info for irrep %i' % iirrep)
                if rstr == '':
                    break

                words = rstr.split()

                self['ncore_dict'][words[0]] = int(words[1])

            self.ostr += 'ncore=%s\n' % str(self['ncore_dict'])


    def parse_ncore_from_cidrtmsin(self, lines):
        irreps = []
        fc_flags = []

        i = 0
        n_irrep = None

        while i < len(lines):
            line = lines[i].strip()
            if "input the number of irreps" in line:
                try:
                    n_irrep = int(line.split()[0])
                except:
                    pass
            if "symmetry labels" in line:
                i += 1
                break
            i += 1

        if n_irrep is None:
            raise ValueError("Could not determine number of irreps")

        for j in range(n_irrep):
            if i >= len(lines):
                break
            token = lines[i].split("/")[0].strip()
            if token:
                irreps.append(token)
            i += 1

        while i < len(lines):
            line = lines[i].strip()
            if "frozen core orbitals" in line:
                clean = line.split("/")[0].strip()
                if clean:
                    fc_flags = list(map(int, clean.split()))
                break
            i += 1

        idx2irrep = {i + 1: ir for i, ir in enumerate(irreps)}

        #print("DEBUG irreps:", irreps)
        #print("DEBUG fc_flags:", fc_flags)

        ncore = {ir: 0 for ir in irreps}

        for i in range(0, len(fc_flags), 2):
            irrep_idx = fc_flags[i]
            ir = idx2irrep.get(irrep_idx)
            if ir is not None:
                ncore[ir] += 1

        return ncore

    def output_options(self):
        if self.ret_yn('Adjust detailed output options?', False):
            self.read_str('Name of the output file', 'output_file', 'summ.txt')
            op = self.ret_str('Output precision, enter as: (<digits>, <dec. digits>)', '(7,3)')
            self.write_option('output_prec', eval(op))
            self.read_str('Format for molden coefficients', 'mcfmt', '% 10E')
            self.read_yn('Print-out sorted by energies?', 'print_sorted', True)

    def get_rassi_list(self):
        print("""
It is assumed that you ran a RASSI job using the TRD1 option
    and copied the TRD2* files to a directory.
Please specify this directory and choose, which files will be analyzed,
    e.g. only transitions to the ground state.
        """)
        ddir = self.ret_str('Directory with the TRD2 files:', 'TRD', True)
        print("%s contains the following files:"%ddir)

        lfiles = sorted(os.listdir(ddir))
        self.print_list(lfiles)
        rstr = self.ret_str('Input indices of required files (separated by spaces)\n')
        dinds = [int(ino)-1 for ino in rstr.replace(',',' ').split()]

        dlist = [os.path.join(ddir,lfiles[dind]) for dind in dinds]
        self.write_option('ana_files', dlist)

def run_theoinp(afile=''):
    wopt = write_options_theo('dens_ana.in')

    if not afile is None:
        wopt.only_at_lists(afile)
        return

    wopt['prop_list'] = []

    wopt.choose_rtype()
    wopt.set_read_options()
    wopt.choose_mo_file()

    if wopt['rtype'] in ['nos']:
        dotden = False
    else:
        dotden = wopt.ret_yn('Analysis of transition density matrices?', True)
        if dotden:
            if wopt.ret_yn('Perform CT number analysis?', True):
                wopt.make_at_lists()
                wopt.set_Om_desc()
                wopt.set_eh_pop()

            wopt.nto_ana()
            wopt.comp_rho0n()
            try:
                if wopt['comp_rho0n'] or wopt['cube_orbitals'] or wopt['comp_p_h_dens']:
                    wopt.read_int('Number of CPUs for orbkit calculations','numproc', 4)
            except error_handler.MsgError:
                pass

            if wopt.ret_yn('Perform exciton analysis?', True):
                wopt.exciton_ana()
            if wopt.ret_yn('Add descriptors for ionic states?', False):
                wopt['prop_list'] += ['QTa', 'QT2', 'LOCa', 'LOC']

    if (not dotden) and (wopt['rtype'] in ['nos', 'colmcscf', 'rassi', 'libwfa']):
        if wopt.ret_yn('Analysis of state density matrices?', True):
            wopt.read_yn('Print out Mulliken populations?', 'pop_ana', True)
            if wopt['rtype'] in ['nos']:
                if wopt.read_yn('Compute number of unpaired electrons?', 'unpaired_ana', True):
                    wopt['prop_list'] += ['nu', 'nunl']
            wopt.AD_ana()
            wopt.BO_ana()
            wopt.ddm_parse()

    wopt.write_prop_list()

    # Program specific input
    if wopt['rtype'] in ['colmrci']:
        wopt.get_ncore()
    elif wopt['rtype'] in ['rassi'] and not wopt['read_libwfa']:
        wopt.get_rassi_list()

    wopt.output_options()

    wopt.flush(lvprt=1, choose_file=True)

    if wopt['rtype'] == 'colmcscf':
        print("\nNow, please run write_den.bash to prepare the MCSCF density matrices!")

class TheodoreInput(Action):

    name = 'theoinp'

    _user_input = """
    # Fragment def (at_lists) using coords from afile
    afile = :: str, optional, alias=a
    """

    _colt_description = "Input generation for TheoDORE"

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..input_options': 'input_options',
            '..lib_struc': 'lib_struc',
            '..error_handler': 'error_handler',
            '..orbkit_interface': 'orbkit_interface',
    })

    def run(afile):
        theo_header.print_header(__class__._colt_description)
        run_theoinp(afile)
