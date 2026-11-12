"""
Analysis routines for transition density matrices.
"""

from __future__ import print_function, division

from . import dens_ana_base, Om_descriptors, lib_mo, error_handler, pop_ana
from . import orbkit_interface, fchk_parser
import numpy
import os

numpy.set_printoptions(precision=6, suppress=True)

class tden_ana(dens_ana_base.dens_ana_base):
    """
    Analysis of transition density matrices.
    """
    # TODO: more efficient treatment for sparse matrices.

#--------------------------------------------------------------------------#
# Print out
#--------------------------------------------------------------------------#

    def print_all_tden(self, lvprt=2):
        """
        Print out all transition density matrices.
        """
        title = "Transition density matrices (MO basis)"
        function = self.print_tden
        self.printer_base(title, function, lvprt)

    def print_tden(self, state, lvprt=2):
        tden = state['tden']
        Om = numpy.dot(tden.flatten(), tden.flatten())

        print("Omega = %10.7f"%Om)
        if lvprt>=2: print(tden)

#---

    def print_all_OmAt(self, lvprt=2):
        """
        Print out Omega matrices with respect to atoms.
        """
        title = "Omega matrices with respect to atoms"
        function = self.print_OmAt
        self.printer_base(title, function, lvprt)

    def print_OmAt(self, state, lvprt=2):
        Om, OmAt = self.ret_Om_OmAt(state)

        print("Omega = %10.7f"%Om)
        if lvprt>=2: print(OmAt)

#---

    def print_all_OmFrag(self, lvprt=2):
        """
        Print out Omega matrices with respect to fragments.
        """
        title = "Omega matrices with respect to fragments"
        function = self.print_OmFrag
        self.printer_base(title, function, lvprt)

    def print_OmFrag(self, state, lvprt=2):
        Om, OmFrag = self.ret_Om_OmFrag(state)

        print("Omega = %10.7f"%Om)
        if lvprt>=2: print(OmFrag)

    def fprint_OmFrag(self, fname="OmFrag.txt"):
        """
        Print a file containing the Omega matrix.
        """
        if self.ioptions['print_sorted']:
            # Sort by excitation energy
            iterlist = sorted(self.state_list, key=lambda state: state['exc_en'])
        else:
            iterlist = self.state_list

        omf = open(fname, 'w')

        omf.write("%i\n"%len(self.ioptions['at_lists']))
        for state in iterlist:
            Om, OmFrag = self.ret_Om_OmFrag(state)
            if Om is None:
                continue
            omf.write("%10s %8.5f"%(state['name'], Om))
            for el in OmFrag.flatten():
                omf.write(" %8.5f"%el)
            omf.write("\n")

        omf.close()
#---

    def print_all_Om_descriptors(self, lvprt=2, desc_list=[]):
        """
        Print Omega descriptors.
        """
        title = "Omega descriptors"
        function = self.print_Om_descriptors
        self.printer_base(title, function, lvprt, desc_list=desc_list)

    def print_Om_descriptors(self, state, lvprt=2, desc_list=[]):
        if not 'Om_desc' in state:
            Om, OmFrag = self.ret_Om_OmFrag(state)
            state['Om_desc'] = Om_descriptors.Om_desc_coll(Om, OmFrag)

        print(self.ret_header_string(desc_list))
        print(state['Om_desc'].ret_val_string(desc_list))

#---

    def print_all_exciton(self, lvprt=2):
        """
        Print out exciton analysis.
        """
        title = "Exciton analysis"
        function = self.print_exciton
        self.printer_base(title, function, lvprt)

    def print_exciton(self, state, lvprt=2):
        Om, OmAt = self.ret_Om_OmAt(state)

        print("RMS e-h sep.: %8.6f Ang"%state['RMSeh'])

#---

    def print_all_eh_pop(self, lvprt=2):
        """
        Print electron/hole population analysis.
        """
        if self.ioptions['eh_pop'] == 0: return

        print("\n*** Electron/hole population analysis ***")

        if 'at_lists' in self.ioptions:
            title = 'Decomposition over fragments'
            function = self.print_eh_Frag
            self.printer_base(title, function, lvprt)
            self.fprint_ehFrag()

        if self.ioptions['eh_pop'] >= 2:
            title = 'Decomposition over individual atoms'
            function = self.print_eh_At
            self.printer_base(title, function, lvprt)

        if self.ioptions['eh_pop'] >= 3:
            title = 'Decomposition over individual basis functions'
            function = self.print_eh_Bas
            self.printer_base(title, function, lvprt)

    def print_eh_Frag(self, state, lvprt=2):
        Om, OmFrag = self.ret_Om_OmFrag(state)

        if OmFrag is None: return

        hpop = numpy.sum(OmFrag, 1)
        epop = numpy.sum(OmFrag, 0)

        for ifrag in range(len(self.ioptions['at_lists'])):
            state["H_%i"%(ifrag+1)] = hpop[ifrag]
            state["E_%i"%(ifrag+1)] = epop[ifrag]

        pop_pr = pop_ana.pop_printer(self.struc)
        pop_pr.add_pop('h+', hpop)
        pop_pr.add_pop('e-', epop)
        pop_pr.add_pop('sum', hpop+epop)
        pop_pr.add_pop('diff', hpop-epop)

        print(pop_pr.ret_table_Frag(self.ioptions['at_lists']))

    def fprint_ehFrag(self, fname="ehFrag.txt"):
        """
        Print a file containing the e/h populations.
        """
        eh_list = []
        for ifrag in range(len(self.ioptions['at_lists'])):
            eh_list.append("H_%i"%(ifrag+1))
            eh_list.append("E_%i"%(ifrag+1))

        ostr = self.ret_summ_table(eh_list)

        open(fname, 'w').write(ostr)
        print("File %s with information about e/h populations written."%fname)

    def print_eh_At(self, state, lvprt=2):
        Om, OmAt = self.ret_Om_OmAt(state)

        hpop = numpy.sum(OmAt, 1)
        epop = numpy.sum(OmAt, 0)

        pop_pr = pop_ana.pop_printer(self.struc)
        pop_pr.add_pop('h+', hpop)
        pop_pr.add_pop('e-', epop)
        pop_pr.add_pop('sum', hpop+epop)
        pop_pr.add_pop('diff', hpop-epop)
        pop_pr.add_pop('trans', state['tpop'])

        print(pop_pr.ret_table())

    def print_eh_Bas(self, state, lvprt=2):
        OmBas = state['OmBas']

        hpop = numpy.sum(OmBas, 1)
        epop = numpy.sum(OmBas, 0)

        pop_pr = pop_ana.pop_printer()
        pop_pr.add_pop('h+', hpop)
        pop_pr.add_pop('e-', epop)
        pop_pr.add_pop('sum', hpop+epop)
        pop_pr.add_pop('diff', hpop-epop)

        print(pop_pr.ret_table())

#--------------------------------------------------------------------------#
# Find data
#--------------------------------------------------------------------------#

    def ret_prop_val(self, prop, state):
        """
        Find the value of a property.
        """
        if prop in state:
            return state[prop]
        else:
            if not 'Om_desc' in state:
                Om, OmFrag = self.ret_Om_OmFrag(state)
                if Om is None: return None

                state['Om_desc'] = Om_descriptors.Om_desc_coll(Om, OmFrag)

            return state['Om_desc'].ret_desc(prop)

#--------------------------------------------------------------------------#
# Operations
#--------------------------------------------------------------------------#

    def compute_all_OmAt(self, fullmat=False):
        """
        Computation of Omega matrices and storage in memory.

        For fullmat==True, also the off-diagonal elements are computed.
        """
        if self.ioptions['print_OmAt'] and os.path.exists('OmAt.npy'):
            Omtmp = numpy.load('OmAt.npy')
            print("Reading Omega matrix from OmAt.npy")
            for i, state in enumerate(self.state_list):
                state['Om']   = numpy.sum(Omtmp[i])
                state['OmAt'] = Omtmp[i]
        else:
            for state in self.state_list:
                Om, OmAt = self.ret_Om_OmAt(state, fullmat)
            if self.ioptions['print_OmAt']:
                numpy.save('OmAt.npy', [state['OmAt'] for state in self.state_list])

        if fullmat:
            self.compute_OmAt_mat()

    def ret_Om_OmAt(self, state, fullmat=False):
        """
        Construction of the Omega matrix with respect to atoms.

        formula=0: Om_mn = (DS)_mn (SD)_mn [JCTC (2012) 8, 2777]
        formula=1: Om_mn = 1/2 (DS)_mn (SD)_mn + 1/2 D_mn (SDS)_mn [JCP (2014), 141, 024106]
        formula=2: Lowdin partitioning
           Om = S^.5 DAO S^.5 = UV^T DMO (UV^T)^T - U, V singular vectors of C
        """
        if 'Om' in state and 'OmAt' in state:
            return state['Om'], state['OmAt']

        formula = self.ioptions.get('Om_formula')

        try:
            D  = state['tden']
        except KeyError:
            return None, None

        print("Computation of Omega matrix ...")

        if formula <= 1:
            OmBas = self.mos.OmBas_Mulliken(D, formula)
            state['tpop'] = 2**.5 * pop_ana.mullpop_ana().ret_pop(D, self.mos)
            state['LOCa'] = 2**.5 * sum(numpy.sqrt(abs(numpy.diag(OmBas))))
        elif formula == 2:
            SDSh = self.mos.lowdin_trans(D)
            if fullmat or self.ioptions['comp_dntos']:
                state['SDSh'] = SDSh
            OmBas = SDSh * SDSh
            state['tpop'] = 2**.5 * pop_ana.pop_ana().ret_pop(D, self.mos, SDSh)
            state['LOCa'] = 2**.5 * numpy.trace(abs(SDSh))
        elif formula == 3: # element-wise multiplication of D and S
            # The results are to be interpreted differently from the above cases.
            # This shows overlap populations of covalent states, and not charge transfer
            # The overlap populations are represented via CT or Om - LOC
            DeS = self.mos.D_element_S(D)
            OmBas = DeS * DeS
            state['tpop'] = 2**.5 * pop_ana.pop_ana().ret_pop(D, self.mos, DeS)
            state['LOCa'] = 2**.5 * numpy.trace(abs(DeS))
        else:
            raise error_handler.MsgError("Om_formula=%i for CT numbers not implemented!"%formula)

        # store OmBas if needed
        if self.ioptions['eh_pop'] >= 3:
            state['OmBas'] = OmBas

        state['LOC'] = numpy.trace(OmBas)
        state['QT2'] = numpy.sum(state['tpop']**2)
        state['QTa'] = numpy.sum(abs(state['tpop']))

        state['Om'] = numpy.sum(OmBas)
        state['OmAt'] = self.mos.comp_OmAt(OmBas)

        return state['Om'], state['OmAt']

    def compute_OmAt_mat(self):
        """
        Compute the full matrix including off-diagonal OmAt elements.
        This is only supported for Lowdin orthogonalization.
        """
        print("Computation of off-diagonal Omega matrices ...")

        if not self.ioptions.get('Om_formula') == 2:
            raise error_handler.MsgError('Only Om_formula==2 supported for full OmAt matrix')

        nstate = len(self.state_list)
        bf_blocks = self.mos.bf_blocks()

        self.Om_mat = numpy.zeros([nstate, nstate], float)
        self.Om_At_mats = [[numpy.zeros([self.mos.num_at, self.mos.num_at]) for i in range(nstate)] for j in range(nstate)]
        for i, istate in enumerate(self.state_list):
            self.Om_mat[i, i] = istate['Om']
            self.Om_At_mats[i][i] = istate['OmAt']
            for j in range(i+1, nstate):
                jstate = self.state_list[j]
                #if ('mult' in istate and 'mult' in jstate) and (istate['mult'] != jstate['mult']): continue
                OmBasIJ = istate['SDSh'] * jstate['SDSh']
                self.Om_mat[i, j] = numpy.sum(OmBasIJ)
                self.Om_mat[j, i] = self.Om_mat[i, j]

                for iat, ist, ien in bf_blocks:
                    for jat, jst, jen in bf_blocks:
                        self.Om_At_mats[i][j][iat, jat] = numpy.sum(OmBasIJ[ist:ien, jst:jen])
                        self.Om_At_mats[j][i][iat, jat] = self.Om_At_mats[i][j][iat, jat]
#---

    def compute_all_OmFrag(self):
        """
        Computation of Omega matrices and storage in memory.
        """
        if 'at_lists' not in self.ioptions:
            print('\n WARNING: at_lists not defined - not computing CT numbers!\n')
            return

        for state in self.state_list:
            Om, OmFrag = self.ret_Om_OmFrag(state, self.ioptions.get('at_lists'))

    def ret_Om_OmFrag(self, state, at_lists=None):
        if at_lists is None:
            try:
                return state['Om'], state['OmFrag']
            except:
                return None, None
                #raise error_handler.MsgError("at_lists not specified")

        Om, OmAt = self.ret_Om_OmAt(state)
        if Om is None:
            return None, None

        state['OmFrag'] = numpy.zeros([len(at_lists), len(at_lists)])

        for A, Aatoms in enumerate(at_lists):
            for B, Batoms in enumerate(at_lists):
                for Aatom in Aatoms:
                    for Batom in Batoms:
                        try:
                            state['OmFrag'][A, B] += state['OmAt'][Aatom-1, Batom-1]
                        except IndexError:
                            print("\n  ERROR: Please check at_lists definition.\n")
                            raise

        return state['Om'], state['OmFrag']

#---
    def compute_all_Phe(self):
        """
        Computation of the deexcitation values <P_he>.
        """
        for i, state in enumerate(self.state_list):
            self.compute_Phe(state)

    def compute_Phe(self, state):
        if 'Phe' in state:
            return

        try:
            D  = state['tden']
        except KeyError:
            return

        if D.shape[0] == D.shape[1]:
            state['Phe'] = numpy.sum(D * D.T) / state['Om']
        elif D.shape[0] < D.shape[1]:
            state['Phe'] = numpy.sum(D[:,:D.shape[0]] * D[:,:D.shape[0]].T) / state['Om']
        else:
            raise error_handler.ElseError('>', 'D.shape')

#---

    def compute_all_NTO(self):
        if len(self.state_list) == 0: return
        if not 'tden' in self.state_list[0]: return

        jmol_orbs = self.ioptions.get('jmol_orbitals')
        if jmol_orbs:
            jmolNTO = lib_mo.jmol_MOs("nto")
            jmolNTO.pre(ofile=self.ioptions.get('mo_file', strict=False))
        cube_ids = []
        for state in self.state_list:
            (U, lam, Vt) = self.ret_NTO(state)
            if jmol_orbs:
                self.export_NTOs_jmol(state, jmolNTO, U, lam, Vt, minlam=self.ioptions['min_occ'])

            if self.ioptions['molden_orbitals']:
                self.export_NTOs_molden(state, U, lam, Vt, minlam=self.ioptions['min_occ'])

            if self.ioptions.get('cube_orbitals'):
                lib_orbkit = orbkit_interface.lib_orbkit(self.ioptions['orbkit_extend'], self.ioptions['orbkit_step'])
                if (lam > self.ioptions['min_occ']).any():
                    cbfid = lib_orbkit.cube_file_creator(state, U, lam, Vt, self.mos,minlam=self.ioptions['min_occ'],numproc=self.ioptions.get('numproc'))
                    cube_ids.append(cbfid)
                else:
                    print(f"WARNING: no singular values above the default min_occ={self.ioptions['min_occ']} found. No NTO cube files will be created for this state.")

        if self.ioptions.get('vmd_ntos'):
            print("VMD network for NTOs")
            lib_orbkit.vmd_network_creator(filename='NTOs',cube_ids=numpy.hstack(cube_ids),isovalue=self.ioptions.get('vmd_ntos_iv'))

        if jmol_orbs:
            jmolNTO.post()

    def compute_p_h_dens(self):
        """
        Computation of electron/hole densities.
        """
        if len(self.state_list) == 0: return
        if not 'tden' in self.state_list[0]: return

        lib_orbkit = orbkit_interface.lib_orbkit(self.ioptions['orbkit_extend'], self.ioptions['orbkit_step'])
        cube_ids = []
        for state in self.state_list:
            (U, lam, Vt) = self.ret_NTO(state)
            if (lam > self.ioptions['min_occ']).any():
                cbfid = lib_orbkit.compute_p_h_dens(state, U, lam, Vt, self.mos, minlam=self.ioptions['min_occ'],numproc=self.ioptions.get('numproc'))
                cube_ids.append(cbfid)
            else:
                print(f"WARNING: no singular values above the default min_occ={self.ioptions['min_occ']} found. No p_h_dens will be created for this state.")

        if self.ioptions.get('vmd_ph_dens'):
            print("VMD network for particle/hole densities")
            lib_orbkit.vmd_network_creator(filename='p_h_dens',cube_ids=numpy.hstack(cube_ids),isovalue=self.ioptions.get('vmd_ph_dens_iv'))

    def compute_rho_0_n(self):
        """
        Computation of transition densities.
        """
        if len(self.state_list) == 0: return
        if not 'tden' in self.state_list[0]: return

        lib_orbkit = orbkit_interface.lib_orbkit(self.ioptions['orbkit_extend'], self.ioptions['orbkit_step'])
        cube_ids = lib_orbkit.compute_rho_0_n(self.state_list,self.mos,numproc=self.ioptions.get('numproc'))
        if self.ioptions.get('vmd_rho0n'):
            print("VMD network for transition densities")
            lib_orbkit.vmd_network_creator(filename='rho0n',cube_ids=numpy.hstack(cube_ids),isovalue=self.ioptions.get('vmd_rho0n_iv'))

    def ret_NTO(self, state):
        if not 'tden' in state: return None, None, None

        # sqrlam contains the squareroot of the singular values lambda as defined in JCP 141, 024106 (2014).
        (U, sqrlam, Vt) = numpy.linalg.svd(state['tden'])
        lam = sqrlam * sqrlam
        lams = lam.sum()

        if lams > 1.e-6:
            state['PRNTO'] = lams * lams / (lam*lam).sum()

        # entanglement entropy
        #   take out the zeros since 0*log(0)=0
        # loglam = numpy.array([0. if lami <= 0. else numpy.log2(lami/lams) for lami in lam])
        # state['S_HE'] = -2.*sum(lam/lams * loglam)
        loglam = numpy.array([0. if lami <= 0. else numpy.log2(lami) for lami in lam])
        state['S_HE'] = -sum(lam * loglam)
        state['Z_HE'] = 2.**(state['S_HE'])

        return U, lam, Vt

    def export_NTOs_jmol(self, state, jmolNTO, U, lam, Vt, mincoeff=0.2, minlam=0.05, pref='NTO', post=''):
        Ut = numpy.transpose(U)
        sname = pref + state['name'].replace('(', '-').replace(')', '-') + post
        jmolNTO.next_set(sname)
        for i, l in enumerate(lam):
            if l < minlam: break

            jmolI = 'mo color blue red\nmo ['
            jmolF = 'mo color orange green\nmo ['

            for occind in (-abs(Ut[i])**2.).argsort():
                occ = Ut[i][occind]
                if abs(occ) < mincoeff: break

                jmolI += ' %.3f %i'%(occ,occind+1)

            for virtind in (-abs(Vt[i])**2.).argsort():
                virt = Vt[i][virtind]
                if abs(virt) < mincoeff: break

                jmolF += ' %.3f %i'%(virt,virtind+1)

            jmolI += ']\n'
            jmolNTO.add_mo(jmolI, "%s_%io"%(sname,i+1), l)
            jmolF += ']\n'
            jmolNTO.add_mo(jmolF, "%s_%iv"%(sname,i+1), l)

    def export_NTOs_molden(self, state, U, lam, Vt, mincoeff=0.2, minlam=0.01, pref='nto', post=''):
        """
        Export the NTOs to a molden file.
        """
        mld_name = '%s_%s%s.mld'%(pref,state['name'].replace('(', '-').replace(')', '-'),post)
        self.mos.export_NTO(lam, U, Vt, mld_name,
                           cfmt=self.ioptions['mcfmt'], occmin=minlam, alphabeta=self.ioptions['alphabeta'])

#---
    def compute_all_DNTO(self):
        """
        Computation of domain NTOs and conditional electron/hole densities.
        """
        if len(self.state_list) == 0: return
        if not 'tden' in self.state_list[0]: return

        print("\n*** Computing domain NTOs ... ***")
        if not self.ioptions.get('Om_formula') == 2:
            raise error_handler.MsgError('Only Om_formula==2 supported for conditional e/h densities.')

        jmol_orbs = self.ioptions['jmol_orbitals']
        if jmol_orbs:
            jmh = lib_mo.jmol_MOs("dnto_hole")
            jmh.pre(ofile=self.ioptions.get('mo_file', strict=False))
            jme = lib_mo.jmol_MOs("dnto_elec")
            jme.pre(ofile=self.ioptions.get('mo_file', strict=False))

        dnto_dens = self.ioptions['comp_dnto_dens']
        if dnto_dens > 0:
            lib_orbkit = orbkit_interface.lib_orbkit(self.ioptions['orbkit_extend'], self.ioptions['orbkit_step'])
            cube_ids = []

        fchk_dens = self.ioptions['fchk_dnto_dens']
        if fchk_dens > 0:
            nb = self.mos.ret_num_bas()
            DNTO_denss = [numpy.zeros([nb, nb], float), numpy.zeros([nb, nb], float)]
            fex = fchk_parser.fchk_export(self.ioptions['rfile'])
        else:
            DNTO_denss = None

        for state in self.state_list:
            print("DNTOs for ", state['name'])

            for A, Aatoms in enumerate(self.ioptions['at_lists']):
                if not self.ioptions['dnto_frags'] == []:
                    if not (A+1) in self.ioptions['dnto_frags']:
                        continue
                #export_opts={'minlam':self.ioptions['min_occ'], 'pref':"DNTO_"}
                export_opts={'minlam':self.ioptions['min_occ']}

                ### conditional hole density ###
                (U, lam, Vt) = self.ret_DNTO_h(state, Aatoms, DNTO_denss)
                export_opts['post'] = "_hole-F%02i"%(A+1)
                if jmol_orbs:
                    self.export_NTOs_jmol(state, jmh, U, lam, Vt, **export_opts)
                if self.ioptions['molden_orbitals']:
                    self.export_NTOs_molden(state, U, lam, Vt, **export_opts)
                if fchk_dens == 1 or fchk_dens == 3:
                    fex.dump_LTmat('%s hole-F%02i Hole Density'%(state['name'], A+1), DNTO_denss[0])
                    fex.dump_LTmat('%s hole-F%02i Electron Density'%(state['name'], A+1), DNTO_denss[1])
                if dnto_dens == 1 or dnto_dens == 3:
                    if sum(lam) < self.ioptions['min_occ']:
                        print("Norm = %.2f < min_occ. Skipping %s ..."%(sum(lam), export_opts['post']))
                    else:
                        N = 1/sum(lam) if self.ioptions['normalize_dnto_dens'] else 1.
                        try:
                            cbfid = lib_orbkit.compute_p_h_dens(state, U, N * lam, Vt,
                            self.mos, numproc=self.ioptions['numproc'], **export_opts)
                            cube_ids.append(cbfid)
                        except:
                            print("... failed.")

                ### conditional electron density ###
                (U, lam, Vt) = self.ret_DNTO_e(state, Aatoms, DNTO_denss)
                export_opts['post'] = "_elec-F%02i"%(A+1)
                if jmol_orbs:
                    self.export_NTOs_jmol(state, jme, U, lam, Vt, **export_opts)
                if self.ioptions['molden_orbitals']:
                    self.export_NTOs_molden(state, U, lam, Vt, **export_opts)
                if fchk_dens >= 2:
                    fex.dump_LTmat('%s elec-F%02i Hole Density'%(state['name'], A+1), DNTO_denss[0])
                    fex.dump_LTmat('%s elec-F%02i Electron Density'%(state['name'], A+1), DNTO_denss[1])
                if dnto_dens >= 2:
                    if sum(lam) < self.ioptions['min_occ']:
                         print("Norm = %.2f < min_occ. Skipping %s ..."%(sum(lam), export_opts['post']))
                    else:
                        N = 1/sum(lam) if self.ioptions['normalize_dnto_dens'] else 1.
                        try:
                            cbfid = lib_orbkit.compute_p_h_dens(state, U, N * lam, Vt,
                            self.mos, numproc=self.ioptions['numproc'], **export_opts)
                            cube_ids.append(cbfid)
                        except:
                            print("... failed.")

        if jmol_orbs:
            jmh.post()
            jme.post()

    def ret_DNTO_h(self, state, Aatoms, DNTO_denss=None):
        # Compute an SVD for the density matrix with hole
        #   coordinates restricted to fragment A
        # TODO: Operate only on the non-zero blocks
        D = state['SDSh']
        DA = numpy.zeros(D.shape, float)
        for iat, ist, ien in self.mos.bf_blocks():
            if iat+1 in Aatoms:
                DA[ist:ien,:] = D[ist:ien,:]

        (U, sqrlam, Vt) = numpy.linalg.svd(self.mos.lowdin_trans(DA, True))
        lam = sqrlam * sqrlam

        if not DNTO_denss is None:
            DNTO_denss[0] = self.mos.lowdin_AO_trans(numpy.dot(DA, DA.T))
            DNTO_denss[1] = self.mos.lowdin_AO_trans(numpy.dot(DA.T, DA))

        return U, lam, Vt

    def ret_DNTO_e(self, state, Aatoms, DNTO_denss=None):
        # Compute an SVD for the density matrix with electron
        #   coordinates restricted to fragment A
        # TODO: Operate only on the non-zero blocks
        D = state['SDSh']
        DA = numpy.zeros(D.shape, float)
        for iat, ist, ien in self.mos.bf_blocks():
            if iat+1 in Aatoms:
                DA[:,ist:ien] = D[:,ist:ien]

        (U, sqrlam, Vt) = numpy.linalg.svd(self.mos.lowdin_trans(DA, True))
        lam = sqrlam * sqrlam

        if not DNTO_denss is None:
            DNTO_denss[0] = self.mos.lowdin_AO_trans(numpy.dot(DA, DA.T))
            DNTO_denss[1] = self.mos.lowdin_AO_trans(numpy.dot(DA.T, DA))

        return U, lam, Vt

    def export_p_h_obs_molden(self, label, om, T, occmin=0.1):
        self.mos.export_MO(om, om, T, label,
           cfmt=self.ioptions['mcfmt'], occmin=self.ioptions['min_occ'])

    #def export_p_h_orbs_jmol(self, state, jmol_ph)

    def compute_rho_0_n(self):
        """
        Computation of transition densities.
        """
        if len(self.state_list) == 0: return
        if not 'tden' in self.state_list[0]: return

        lib_orbkit = orbkit_interface.lib_orbkit(self.ioptions['orbkit_extend'], self.ioptions['orbkit_step'])
        cube_ids = lib_orbkit.compute_rho_0_n(self.state_list,self.mos,numproc=self.ioptions.get('numproc'))
        if self.ioptions.get('vmd_rho0n'):
            print("VMD network for transition densities")
            lib_orbkit.vmd_network_creator(filename='rho0n',cube_ids=numpy.hstack(cube_ids),isovalue=self.ioptions.get('vmd_rho0n_iv'))

#---

    def analyze_excitons(self, exciton_ana):
        for state in self.state_list:
            Om, OmAt = self.ret_Om_OmAt(state)
            if Om is None: continue

            state['RMSeh'] = exciton_ana.ret_RMSeh(Om, OmAt)
            state['MAeh']  = exciton_ana.ret_MAeh(Om, OmAt)
            state['K2']    = exciton_ana.ret_K2(OmAt)

            try:
                state['rTD']  = exciton_ana.ret_rTD(state['QT2'], state['tpop'])
                state['Vint'] = exciton_ana.ret_Vinter(state['QT2'], state['tpop'])
                state['Vdiag'] = exciton_ana.ret_Vdiag(state['QT2'], state['tpop'])
                state['J2'] = exciton_ana.ret_J2(state['QT2'], state['tpop'])
            except KeyError:
                pass
#---

    def compute_es2es_tden(self, iref=1):
        """
        Compute approximate transition density with respect to a reference excited state. This overwrites 'tden' and 'exc_en'
        This is computed analogously to the electron/hole densities:
        D^IJ = (D^I0)^T * D^J0 - D^J0 * (D^I0)^T
        """
        refstate = self.state_list[iref-1]
        tdenI = refstate['tden']
        enI   = refstate['exc_en']
        print("Computing approximate state-to-state transition densities")
        print(" ... reference state: %s \n"%refstate['name'])

        for state in self.state_list:
            tdenJ = state['tden']

            DIJ_elec = numpy.dot(tdenI.T, tdenJ)
            DIJ_hole = numpy.dot(tdenJ, tdenI.T)

            state['tden'] = DIJ_elec
            state['tden'][:DIJ_hole.shape[0], :DIJ_hole.shape[1]] -= DIJ_hole
            state['tden'] /= 2**.5

            state['exc_en'] = state['exc_en'] - enI
            del state['osc_str']

        del self.state_list[iref-1]
