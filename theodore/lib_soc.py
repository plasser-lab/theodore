"""
Routines related to spin-orbit coupling are collected here.
"""

from __future__ import print_function, division

from . import lib_tden, file_parser, lib_mo, lib_struc, units, error_handler, lib_exciton
import numpy
import copy

class tden_ana_soc(lib_tden.tden_ana):
    """
    Analysis of transition density matrices for spin-orbit coupled states.
    """

    def read_dens(self):
        """
        Read the (transition) density matrices and some supplementary information.
        """
        rtype = self.ioptions.get('rtype').lower()

        if rtype == 'adf':
                self.mos = lib_mo.MO_set_adf(file=self.ioptions.get('rfile'))
                self.mos.read(lvprt=1)
                self.read2_mos()
                self.struc = lib_struc.structure()
                self.struc.read_at_dicts(self.mos.at_dicts)

                self.state_list_mch, self.state_list_soc = file_parser_adf_soc(self.ioptions).read_soc(self.mos)
                self.state_list = self.state_list_mch

                # deactivate print out that is not possible because of the use of STOs
                self.ioptions['jmol_orbitals'] = False
                self.ioptions['molden_orbitals'] = False
        else:
            raise error_handler.ElseError(rtype, 'rtype for tden_ana_soc')

    def soc_transform(self):
        """
        Transform the data to the spin-orbit coupled representation
        """
        print("Transforming the Omega matrices to the diagonal representation")
        nsing = len(self.state_list_soc[1]['coeffS'])
        nat = self.mos.num_at
        s2  = numpy.sqrt(.5)
        for state in self.state_list_soc[1:]:
            # add up the Omega matrices into the correct spin blocks
            # alpha-alpha
            OmAt_tmp = numpy.zeros([nat,nat], complex)
            coeffs = (s2*state['coeffS']).tolist() + (s2*state['coeffT'][:,0]).tolist()
            for i, ic in enumerate(coeffs):
                for j, jc in enumerate(coeffs):
                    OmAt_tmp  += ic.conjugate() * jc * self.Om_At_mats[i][j]
            state['OmAt_aa']  = numpy.real(OmAt_tmp)
            state['Om_aa']    = numpy.sum(state['OmAt_aa'])
            assert numpy.sum(numpy.imag(OmAt_tmp)**2) < 1.E-12

            # beta-beta
            OmAt_tmp = numpy.zeros([nat,nat], complex)
            coeffs = (s2*state['coeffS']).tolist() + (-s2*state['coeffT'][:,0]).tolist()
            for i, ic in enumerate(coeffs):
                for j, jc in enumerate(coeffs):
                    OmAt_tmp += ic.conjugate() * jc * self.Om_At_mats[i][j]
            state['OmAt_bb'] = numpy.real(OmAt_tmp)
            state['Om_bb']    = numpy.sum(state['OmAt_bb'])
            assert numpy.sum(numpy.imag(OmAt_tmp)**2) < 1.E-12

            # alpha-beta
            OmAt_tmp = numpy.zeros([nat,nat], complex)
            coeffs = [0. for i in range(nsing)] + state['coeffT'][:,1].tolist()
            for i, ic in enumerate(coeffs):
                for j, jc in enumerate(coeffs):
                    OmAt_tmp += ic.conjugate() * jc * self.Om_At_mats[i][j]
            state['OmAt_ab'] = numpy.real(OmAt_tmp)
            state['Om_ab']    = numpy.sum(state['OmAt_ab'])
            assert numpy.sum(numpy.imag(OmAt_tmp)**2) < 1.E-12

            # beta-alpha
            OmAt_tmp = numpy.zeros([nat,nat], complex)
            coeffs = [0. for i in range(nsing)] + state['coeffT'][:,2].tolist()
            for i, ic in enumerate(coeffs):
                for j, jc in enumerate(coeffs):
                    OmAt_tmp += ic.conjugate() * jc * self.Om_At_mats[i][j]
            state['OmAt_ba'] = numpy.real(OmAt_tmp)
            state['Om_ba']    = numpy.sum(state['OmAt_ba'])
            assert numpy.sum(numpy.imag(OmAt_tmp)**2) < 1.E-12

            state['Om'] = state['Om_aa'] + state['Om_bb'] + state['Om_ab'] + state['Om_ba']
            state['OmAt'] = state['OmAt_aa'] + state['OmAt_bb'] + state['OmAt_ab'] + state['OmAt_ba']

    def print_info(self, pre):
        """
        Print info for different types of states.
        """
        if pre == 'mch':
            self.state_list = self.state_list_mch
            header = "  Analysis of original (MCH) states"
        elif pre == 'aa':
            self.state_list = copy.deepcopy(self.state_list_soc)
            header = "  Analysis of alpha-alpha components"
            for state in self.state_list[1:]:
                state['Om']   = state['Om_aa']
                state['OmAt'] = state['OmAt_aa']
        elif pre == 'bb':
            self.state_list = copy.deepcopy(self.state_list_soc)
            header = "  Analysis of beta-beta components"
            for state in self.state_list[1:]:
                state['Om']   = state['Om_bb']
                state['OmAt'] = state['OmAt_bb']
        elif pre == 'ab':
            self.state_list = copy.deepcopy(self.state_list_soc)
            header = "  Analysis of alpha-beta components"
            for state in self.state_list[1:]:
                state['Om']   = state['Om_ab']
                state['OmAt'] = state['OmAt_ab']
        elif pre == 'ba':
            self.state_list = copy.deepcopy(self.state_list_soc)
            header = "  Analysis of beta-alpha components"
            for state in self.state_list[1:]:
                state['Om']   = state['Om_ba']
                state['OmAt'] = state['OmAt_ba']
        elif pre == 'soc':
            self.state_list = self.state_list_soc
            header = "  Analysis of spin-orbit coupled (diagonal) states"
            if not 'Om_aa' in self.ioptions['prop_list']:
                self.ioptions['prop_list'] += ['Om_aa', 'Om_ab']
        else:
            raise error_handler.ElseError(pre, 'pre')

        print()
        print("=====================================================")
        print(header)
        print("=====================================================")

        self.ioptions['output_file'] = '%s_tden_summ.txt'%pre
        if 'at_lists' in self.ioptions:
            self.compute_all_OmFrag()
            if self.ioptions['print_OmFrag']:
                self.fprint_OmFrag('%s_OmFrag.txt'%pre)
        if 'RMSeh' in self.ioptions.get('prop_list') or 'MAeh' in self.ioptions.get('prop_list') or 'K2' in self.ioptions.get('prop_list'):
            exca = lib_exciton.exciton_analysis()
            exca.get_distance_matrix(self.struc)
            self.analyze_excitons(exca)
        self.print_all_eh_pop()
        self.print_summary()

class file_parser_adf_soc(file_parser.file_parser_adf):
    """
    Read ADF TDDFT SOC information from the TAPE21 file.
    """
    def read_soc(self, mos):
        try:
            from scm.plams import KFFile
        except ImportError:
            from kf import kffile as KFFile
        rfile = KFFile(self.ioptions['rfile'])

        state_list = self.read(mos)

        print("Reading spin-orbit information from ADF...")
        excs = rfile.read('Excitations SO A', 'excenergies') * units.energy['eV']
        oscs = rfile.read('Excitations SO A', 'oscillator strengths')
        nstate = len(excs)

        state_list_soc = [{} for istate in range(nstate)]
        for istate in range(nstate):
            state = state_list_soc[istate]
            state['state_ind'] = istate
            state['exc_en'] = excs[istate]
            state['osc_str'] = oscs[istate]
            state['name'] = 'E%i'%(state['state_ind'])

        real_tri=rfile.read('Excitations SO A','SOmat-R')
        imag_tri=rfile.read('Excitations SO A','SOmat-I')
        Hsoc = numpy.zeros([nstate, nstate], complex)
        x=0
        y=0
        for i in range(len(real_tri)):
            if abs(real_tri[i])<1e-15: real_tri[i]=0.
            if abs(imag_tri[i])<1e-15: imag_tri[i]=0.
            Hsoc[x, y]=real_tri[i] + (0+1j)*imag_tri[i]
            Hsoc[y, x]=real_tri[i] + (0-1j)*imag_tri[i]
            x+=1
            if x>y:
                y+=1
                x=0

        excs_tmp, U = numpy.linalg.eigh(Hsoc)
        assert max(excs_tmp * units.energy['eV'] - excs)<1.E-12, 'Inconsistent Hsoc'

        numpy.savetxt('U.txt', (numpy.conj(U)*U).real, fmt='%.4f', delimiter=' ')

        print("Ground state: <S0|E0>^2 = %.5f"%abs(U[-1,0])**2)

        for a, astate in enumerate(state_list_soc[1:]):
            astate['coeffS'] = numpy.zeros(self.nsing, complex)
            astate['coeffT'] = numpy.zeros([self.ntrip, 3], complex)

            I = 0
            for i in range(self.nsing):
                astate['coeffS'][i] = U[I,a+1]
                I+=1
            for i in range(self.ntrip):
                astate['coeffT'][i,:] = U[I:I+3,a+1]
                I+=3

        rfile.close()

        return state_list, state_list_soc
