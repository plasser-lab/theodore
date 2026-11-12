from __future__ import print_function, division
from . import file_parser, lib_mo, error_handler, cclib_interface, fchk_parser, units, lib_struc
import numpy

class dens_ana_base:
    """
    Base class for density matrix analysis.
    """
    def __init__(self, ioptions):

        # state_list contains all the information about the states:
        #    quantities that are parsed from files as well as computed quantities
        self.state_list = []

        self.ioptions = ioptions

#--------------------------------------------------------------------------#
# Input
#--------------------------------------------------------------------------#

    def read_mos(self, lvprt=1, spin=0):
        """
        Read MOs from a separate file, which is given in Molden format.
        """
        rtype = self.ioptions.get('rtype')
        if rtype=='tddftb':
           self.mos = lib_mo.MO_set_tddftb(file=self.ioptions.get('mo_file'))
        else:
           self.mos = lib_mo.MO_set_molden(file=self.ioptions.get('mo_file'))
        self.mos.read(lvprt=lvprt, spin=spin)
        self.read2_mos(lvprt)

    def read2_mos(self, lvprt=1):
        self.mos.comp_inv_lowdin(self.ioptions['Om_formula'], lvprt)
        self.num_mo  = self.mos.ret_num_mo()
        self.num_bas = self.mos.ret_num_bas()

    def read_dens(self):
        """
        Read the (transition) density matrices and some supplementary information.
        """
        rtype = self.ioptions.get('rtype').lower()
        if self.ioptions['read_libwfa']: self.mos = None

        if rtype=='ricc2':
            self.state_list = file_parser.file_parser_ricc2(self.ioptions).read(self.mos)
        elif rtype in ['tddft', 'escf', 'tmtddft']:
            self.state_list = file_parser.file_parser_escf(self.ioptions).read(self.mos)
        elif rtype=='libwfa':
            self.state_list = file_parser.file_parser_libwfa(self.ioptions).read()
        elif rtype=='qcadc':
            self.state_list = file_parser.file_parser_qcadc(self.ioptions).read()
        elif rtype=='qctddft':
            self.state_list = file_parser.file_parser_qctddft(self.ioptions).read(self.mos)
        elif rtype=='fchk':
            self.mos = fchk_parser.MO_set_fchk(file=self.ioptions.get('rfile'), read=True)
            self.read2_mos()
            self.state_list = fchk_parser.file_parser_fchk(self.ioptions).read(self.mos)
            self.mos.write_molden_file(fname='MOs.mld')
            self.ioptions['mo_file'] = 'MOs.mld'
        elif rtype in ['mcscf', 'colmcscf']:
            self.state_list = file_parser.file_parser_col_mcscf(self.ioptions).read(self.mos)
        elif rtype in ['mrci', 'colmrci']:
            self.state_list = file_parser.file_parser_col_mrci(self.ioptions).read(self.mos)
        elif rtype in ['rassi', 'molcas']:
            self.state_list = file_parser.file_parser_rassi(self.ioptions).read(self.mos)
        elif rtype == 'terachem':
            self.state_list = file_parser.file_parser_terachem(self.ioptions).read(self.mos)
        elif rtype == 'adf':
            self.mos = lib_mo.MO_set_adf(file=self.ioptions.get('rfile'))
            self.mos.read(lvprt=1)
            self.read2_mos()
            self.state_list = file_parser.file_parser_adf(self.ioptions).read(self.mos)

            # deactivate print out that is not possible because of the use of STOs
            self.ioptions['jmol_orbitals'] = False
            self.ioptions['molden_orbitals'] = False
        elif rtype=='tddftb':
            self.state_list = file_parser.file_parser_tddftb(self.ioptions).read(self.mos)
        elif rtype=='dftmrci':
            self.state_list = file_parser.file_parser_dftmrci(self.ioptions).read(self.mos)
        elif rtype=='onetep':
            self.mos = lib_mo.MO_set_onetep(file=self.ioptions.get('rfile'))
            self.mos.read(lvprt=1)
            self.num_mo  = self.mos.ret_num_mo()
            self.num_bas = self.mos.ret_num_bas()
            self.state_list = file_parser.file_parser_onetep(self.ioptions).read(self.mos)
        elif rtype == 'nos':
            self.ioptions['NO_ana'] = False # Do not do explicit NO analysis
            self.state_list = file_parser.file_parser_nos(self.ioptions).read(self.mos)
        elif rtype == 'orca':
            self.state_list = file_parser.file_parser_orca(self.ioptions).read(self.mos)
        elif rtype == 'orcacas':
            self.state_list = file_parser.file_parser_orca_CAS(self.ioptions).read(self.mos)
        elif rtype in ['cclib', 'gamess']:
            # these are parsed with the external cclib library
            ccli = cclib_interface.file_parser_cclib(self.ioptions)

            if not hasattr(self, 'mos'):
                errcode = ccli.check(maxerr=1)
                print()
                self.mos = ccli.read_mos()

            self.read2_mos()
            self.state_list = ccli.read(self.mos)
            self.struc = ccli.ret_struc()
        else:
            raise error_handler.ElseError(rtype, 'rtype')

        self.extra_info()
        if not self.ioptions['ana_states'] == []:
            self.state_list = self.select_states(self.ioptions['ana_states'], self.state_list)

    def extra_info(self, lvprt=1):
        for state in self.state_list:
            try:
                state['lam'] = units.energy['nm'] / (state['exc_en'] / units.energy['eV'])
                state['rcm']  = state['exc_en'] / units.energy['eV'] * units.energy['rcm'] / 1000
            except ZeroDivisionError:
                pass

        if 'coor_file' in self.ioptions:
            if lvprt >= 1:
                print(("\nReading structure from coor_file %s"%self.ioptions['coor_file']))
            self.struc = lib_struc.structure()
            self.struc.read_file(self.ioptions['coor_file'], self.ioptions['coor_format'])
        elif self.ioptions['rtype'].lower() in ['cclib', 'gamess']:
            if lvprt >= 1:
                print("\nUsing cclib structure")
        elif self.ioptions['rtype'].lower() in ['adf']:
            if lvprt >= 1:
                print("\nUsing ADF structure")
            self.struc = lib_struc.structure()
            self.struc.read_at_dicts(self.mos.at_dicts)
        elif 'mo_file' in self.ioptions:
            if lvprt >= 1:
                print(("\nReading structure from mo_file %s"%self.ioptions['mo_file']))
            self.struc = lib_struc.structure()
            self.struc.read_at_dicts(self.mos.at_dicts)
        else:
            self.struc = None

        if not self.struc == None:
            nquad  = self.struc.ret_nuc_multipole(2)
            nquadt = nquad[0] + nquad[1] + nquad[2]
            for state in self.state_list:
                # Quadrupole moment in Buckingham = Debye-Ang
                if 'r2x' in state:
                    # Traceless
                    # state['Qxx'] = (3 * (-state['r2x'] + nquad[0]) + state['r2'] - nquadt) * units.dipole['D'] * units.length['A']
                    # state['Qyy'] = (3 * (-state['r2y'] + nquad[1]) + state['r2'] - nquadt) * units.dipole['D'] * units.length['A']
                    # state['Qzz'] = (3 * (-state['r2z'] + nquad[2]) + state['r2'] - nquadt) * units.dipole['D'] * units.length['A']
                    # Not traceless
                   state['Qxx'] = (-state['r2x'] + nquad[0]) * units.dipole['D'] * units.length['A']
                   state['Qyy'] = (-state['r2y'] + nquad[1]) * units.dipole['D'] * units.length['A']
                   state['Qzz'] = (-state['r2z'] + nquad[2]) * units.dipole['D'] * units.length['A']

            if lvprt >= 1:
                num_at = self.struc.ret_num_at()
                print("Number of atoms: %i"%num_at)
                print("Composition: %s\n"%self.struc.ret_at_list_composition(list(range(1, num_at+1))))

    def select_states(self, ana_states, state_list):
        """
        Analyze only the states given in 'ana_states'
        """
        ret_list = []
        for anas in ana_states:
            ret_list.append(state_list[anas-1])

        return ret_list

#--------------------------------------------------------------------------#
# Output
#--------------------------------------------------------------------------#

    def printer_base(self, title, function, lvprt=2, **kwargs):
        """
        General print-out of properties.
        """
        print()
        print(title)

        for state in self.state_list:
            #print '%i%s'%(state['state_ind'],state['irrep'])
            print(state['name'])
            function(state, lvprt, **kwargs)

#---

    def print_summary(self):
        """
        Print a summary with information specified in prop_list.
        """
        ostr = self.ret_summ_table(self.ioptions.get('prop_list'))

        print("\n" + ostr)

        if 'output_file' in self.ioptions:
            ofile = self.ioptions.get('output_file')
            print("Final output copied to %s"%ofile)
            open(ofile, 'w').write(ostr)

    def ret_summ_table(self, prop_list):
        width, ndec = self.ioptions.get('output_prec')
        oformat = '%% %i.%if'%(width, ndec)

        hstr  = '%-10s'%'state' + '%*s'%(width+1, 'dE(eV)') + '%*s'%(width-2, 'f')
        hstr +=self.ret_header_string(prop_list, width)

        prt_list = []
        for state in self.state_list:
            vstr  = '%-10s'%state['name'][-10:]
            vstr += oformat%state['exc_en']
            try:
                vstr += (oformat%state['osc_str']).replace(oformat%(-1.),'%*s'%(width, '-'))
            except KeyError:
                vstr += '%*s'%(width, '-')

            vstr += self.ret_val_string(prop_list, state, oformat)

            prt_list.append([state['exc_en'], vstr])

        if self.ioptions['print_sorted']: prt_list.sort()

        ostr  = hstr + "\n"
        ostr += len(hstr) * '-' + "-\n"
        for en, vstr in prt_list:
            ostr += vstr + "\n"
        return ostr

    def ret_header_string(self, prop_list, width=7):
        ret_str = ''

        for prop in prop_list:
            ret_str += '%*s'%(width, prop)

        return ret_str

    def ret_val_string(self, prop_list, state, oformat='% 7.3f'):
        ret_str = ''

        for prop in prop_list:
            val = self.ret_prop_val(prop, state)

            if val == None:
                ret_str += (len(oformat%(0.0))-1)*' ' + '-'
            else:
                try:
                    ret_str += oformat%val
                except:
                    ret_str += ' ' + val

        return ret_str

    def ret_prop_val(self, prop, state):
        """
        Find the value of a property.
        -> This is overloaded in lib_tden.py
        """
        if prop in state:
            return state[prop]
        else:
            return None
