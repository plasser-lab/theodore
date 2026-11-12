"""
Parsing of files produced by different quantum chemical programs.
"""

from __future__ import print_function, division
from . import units, lib_mo, error_handler
import numpy
import os, struct


class file_parser_base:
    def __init__(self, ioptions):
        self.ioptions = ioptions

    def read(self):
        state_list = []

        return state_list

    def init_den(self, mos, rect=False):
        """
        Initialize an empty (transition) density matrix.
        rect=True specifies that only the the occ x (occ + virt) block is explicitly constructed.
           This allows to treat CIS like theories efficiently without many other changes in the code.
        """
        num_mo = mos.ret_num_mo()

        if not rect:
            return numpy.zeros([num_mo, num_mo])
        else:
            nocc  = mos.ret_ihomo() + 1
            nvirt = num_mo - nocc
            return numpy.zeros([nocc, num_mo])

    def dens_stat(self, dens, lvprt=1):
        """
        Print statistics about density matrix.
        """
        if lvprt <= 0: return

        print(" Density matrix statistics")
        print(" Trace: % .6f"%numpy.trace(dens))

        print(" Largest elements:")
        for ind in reversed(numpy.argsort((dens*dens).flatten())):
            i = ind // dens.shape[1]
            j = ind %  dens.shape[1]
            print("  (%4i/%4i): % .6f"%(i, j, dens[i, j]))

            if dens[i, j] * dens[i, j] < 0.1:
                break

        if lvprt >= 3:
            print(dens)

    def parse_key(self, state, key, line, search_string, ind=-1, rfile=None, rmatrix=False, not_string=None):
        """
        Find search_string in the specified line and set it as state[key].
        If rfile is given, then the next line is searched for x,y,z components.
        If rmatrix==True, then a tensor is read from the next three lines.
        """
        translate_str = ',[]|'
        if search_string in line:
            if not not_string is None:
                if not_string in line: return

            words = self.delete_chars(line.strip(search_string), translate_str)
            try:
                state[key] = float(words[ind])
            except:
                print(' WARNING: Error in parse_key for %s of state %s'%(key, state['name']))
                print(line)

            if not rfile is None:
                if rmatrix:
                    line = next(rfile)
                    try:
                        words=self.delete_chars(next(rfile), translate_str)
                        state['%sxx'%key] = float(words[-3])
                        state['%sxy'%key] = float(words[-2])
                        state['%sxz'%key] = float(words[-1])
                        words=self.delete_chars(next(rfile), translate_str)
                        state['%syx'%key] = float(words[-3])
                        state['%syy'%key] = float(words[-2])
                        state['%syz'%key] = float(words[-1])
                        words=self.delete_chars(next(rfile), translate_str)
                        state['%szx'%key] = float(words[-3])
                        state['%szy'%key] = float(words[-2])
                        state['%szz'%key] = float(words[-1])
                    except:
                        print(' Cannot parse matrix for %s of state %s'%(key, state['name']))
                else:
                    words=self.delete_chars(next(rfile), translate_str)
                    try:
                        state['%sx'%key] = float(words[-3])
                        state['%sy'%key] = float(words[-2])
                        state['%sz'%key] = float(words[-1])
                    except:
                        print(' WARNING: Error in parse_key for %s of state %s'%(key, state['name']))
                        print(line)
                        print(words)

    def delete_chars(self, line, delete, lmap=None):
        """
        Delete the characters in <delete> from line, split, and optionally map.
        """
        retline = line
        for char in delete:
            retline = retline.replace(char, '')

        return retline.split()

    def sym_split(self, sym):
        """
        Split an MO "sym" label into the index and irrep.
        "51b1u" -> [51, "b1u"]
        """
        st = 1000
        if sym.find('a') != -1: st = min(st, sym.find('a'))
        if sym.find('b') != -1: st = min(st, sym.find('b'))
        if sym.find('e') != -1: st = min(st, sym.find('e'))
        if sym.find('t') != -1: st = min(st, sym.find('t'))

        if (st==1000): raise error_handler.MsgError("sym_split: %s"%sym)

        return (int(sym[:st]), sym[st:])

class configuration:
    def __init__(self):
        self.occ    = "1x" # label of the occupied (hole) orbital
        self.virt   = "1y" # label of the virtual (particle) orbital
        self.coeff  = 0.
        self.weight = 0.

    def read_ricc2_line(self, line):
        try:
            cols = line.split('|')

            words = cols[1].split()
            self.occ = words[0] + words[1]

            words = cols[2].split()
            self.virt = words[0] + words[1]

            words = cols[3].split()
            self.coeff  = float(words[0])
            self.weight = float(words[1])
        except:
            print("\n ERROR when parsing line:")
            print(line)
            raise

#--------------------------------------------------------------------------#
# Implementations
#--------------------------------------------------------------------------#

class file_parser_ricc2(file_parser_base):
    """
    Turbomole ricc2
    """
    def read(self, mos):
        state_list = self.ret_conf_ricc2(rfile=self.ioptions.get('rfile'))

        ihomo  = mos.ret_ihomo()

        for state in state_list:
            state['name'] = '%i(%s)%s'%(state['state_ind'], state['mult'], state['irrep'])
            state['tden'] = self.init_den(mos, rect=True)

            if self.ioptions.get('read_binary'):
                self.set_tden_bin(state, mos)
            else:
                self.set_tden_conf(state, mos)

        if self.ioptions.get('read_binary'):
               mos.symsort(self.ioptions['irrep_labels'], self.ioptions['Om_formula'])
               if self.ioptions['jmol_orbitals']:
                    print(" \nWARNING: jmol_orbitals not possible with read_binary. Use molden_orbitals instead!")
                    self.ioptions['jmol_orbitals'] = False

        return state_list

    def set_tden_conf(self, state, mos):
        """
        Set the transition density matrix elements acoording to the parsed configurations.
        """
        for conf in state['char']:
            # look for the orbital in the MO file. This makes sure that the procedure works even when the MOs are reordered because of symmetry.
            iocc  = mos.syms.index(conf.occ)

            #print conf.virt, mos.syms.index(conf.virt), mos.syms.index(conf.virt) - ihomo, conf.coeff
            ivirt = mos.syms.index(conf.virt)

            try:
                state['tden'][iocc,ivirt] = conf.coeff
            except:
                print("\nERROR setting %s->%s (%i->%i) transition"%(conf.occ, conf.virt, iocc+1, ivirt+1))
                print("Did you delete the line")
                print("       implicit core=   x virt=    x")
                print("  from the control file before running tm2molden?\n")
                raise

    def set_tden_bin(self, state, mos, lvprt=1):
        """
        Set the transition density matrix elements acoording to the binary files with state information.
        """
        iirrep=self.ioptions['irrep_labels'].index(state['irrep']) + 1
        smult=state['mult']
        istate=state['state_ind']
        CCfilen = 'CCRE0{0:->2d}{1:->3s}{2:->4d}'.format(iirrep, smult, istate)

        if lvprt >= 1: print('Reading binary file %s ...'%CCfilen)

        CCfile = open(CCfilen, 'rb')
        CCfile.read(8)
        method = struct.unpack('8s', CCfile.read(8))[0]
        CCfile.read(8)
        nentry = struct.unpack('i', CCfile.read(4))[0]
        CCfile.read(4)

        # The format was changed for TM Version >= 7.2
        vcheck = struct.unpack('l', CCfile.read(8))[0]
        if vcheck == 0:
            if lvprt >= 1:
                print("  Assuming TM Version >= 7.2")
            CCfile.read(16)
        else:
            if lvprt >= 1:
                print("  Assuming TM Version <= 7.2")
            CCfile.read(8)

        num_mo = mos.ret_num_mo()
        nocc  = mos.ret_ihomo() + 1
        nvirt = num_mo - nocc

        assert(nentry % nvirt == 0)
        nact = nentry // nvirt
        nfrzc = nocc - nact

        if lvprt >= 1:
            print('  Method: %s, number of entries: %i'%(method, nentry))
            print('  num_mo: %i, nocc: %i, nvirt: %i, nact: %i, nfrzc: %i'%(num_mo, nocc, nvirt, nact, nfrzc))
        if lvprt >= 2:
            print('  Syms:', mos.syms)

        if nfrzc > 0:
            raise error_handler.MsgError("""Frozen core orbitals detected in the Molden file!
In the case of read_binary=True, do not delete the line
       implicit core=   x virt=    x
from the control file.""")

        # write the collected data into the correct block of the 1TDM
        coeff = struct.unpack(nentry*'d', CCfile.read(nentry*8))
        state['tden'][nfrzc:nocc, nocc:num_mo] = numpy.reshape(coeff, [nact, nvirt])

        lbytes = struct.unpack('4s', CCfile.read(4))
        if lvprt >= 2:
            print('   Last four bytes:', lbytes)

        if not CCfile.read(1) == b'':
            raise error_handler.MsgError('parsing file %s'%CCfilen)

        if lvprt >= 3:
            print('parsed tden:')
            print(state['tden'])

    def ret_conf_ricc2(self, rfile='ricc2.out'):
        """
        Return information about configurations in a Turbomole calculation.
        """
        ret_list = []
        curr_osc = 0 # running index for oscillator strengths
        section  = '' # which section of the file
        lines = open(rfile,'r')
        while True: # loop over all lines
          try:
            line = next(lines)
          except StopIteration:
            print("Finished parsing file %s"%rfile)
            break
          else:
            if (' sym | multi | state' in line) and ('excitation energies' in line):
                # parse the multiplicity from here
                multis = []
                line = next(lines)
                line = next(lines)
                line = next(lines)
                while True:
                    line = next(lines)
                    if '---' in line: continue
                    if '===' in line: break

                    words = self.delete_chars(line, ['|'])

                    multis.append(words[1])

            elif 'Energy:   ' in line:
                ret_list.append({})
                words = line.split()
                ret_list[-1]['exc_en'] = float(words[3])

                for i in range(3): line = next(lines)

                words = line.split() # type: RE0
                ret_list[-1]['irrep'] = words[4]
                ret_list[-1]['state_ind'] = int(words[6])
                ret_list[-1]['mult'] = multis.pop(0)

                for i in range(3): line = next(lines)

                ret_list[-1]['char'] = []
                while True: # loop over information about this excitation
                    line = next(lines)

                    if '===' in line: break

                    conf = configuration()
                    conf.read_ricc2_line(line)
                    ret_list[-1]['char'].append(conf)

            elif 'oscillator strength (length gauge)' in line:
                words = line.split()
                if curr_osc < len(ret_list):
                    ret_list[curr_osc]['osc_str'] = float(words[5])
                elif curr_osc == len(ret_list):
                    print("\nWARNING: More oscillator strengths than singlet states!")
                    print("   This is ok if you are computing state-to-state properties.\n")
                curr_osc+=1

            elif curr_osc < len(ret_list):
                state = ret_list[curr_osc]
                self.parse_key(state, 'Tmux', line, 'xdiplen', 1, not_string='rel')
                self.parse_key(state, 'Tmuy', line, 'ydiplen', 1, not_string='rel')
                self.parse_key(state, 'Tmuz', line, 'zdiplen', 1, not_string='rel')

            elif 'real representations' in line:
                self.ioptions['irrep_labels'] = line.split()[6:]

            elif 'COSMO-ADC(2) energy differences' in line:
                print('\nReading COSMO-ADC(2) energies ...')
                for i in range(4): line = next(lines)

                istate = 0
                while True: # loop over information about this excitation
                    line = next(lines)

                    if '---' in line: continue
                    if '===' in line: break

                    words = line.replace('|','').split()
                    irrep, mult, state_ind, exc_en = words[0].replace('*',''), words[1], int(words[2]), float(words[-1])
                    state = ret_list[istate]
                    if state['irrep'] == irrep and state['mult'] == mult and state['state_ind'] == state_ind:
                        state['exc_en'] = exc_en
                    else:
                        print(('Data not matching:', state['irrep'], irrep, state['mult'], mult, state['state_ind'], state_ind))
                        raise error_handler.MsgError('COSMO-ADC: data not matching')

                    istate += 1

        return ret_list

#---

class file_parser_escf(file_parser_base):
    """
    Turbomole TDDFT.

    Updated by M.Romanelli (May 2026):  unrestricted calculations are supported now.
    """
    def read(self, mos):
        state_list, restr, act_occ, act_virt = self.ret_conf_tddft(rfile=self.ioptions.get('rfile'))

        #print("DEBUG:",restr,act_virt,act_occ)
        #print("DEBUG:",state_list)

        if restr:  #RKS/RHF response calculations
            nocc={}
            nvirt={}
            for state in state_list:
                state['name'] = '%i%s'%(state['state_ind'],state['irrep'])
                state['tden'] = self.init_den(mos, rect=True)
                occmap  = []
                virtmap = []
                for iorb, sym in enumerate(mos.syms):
                    if state['irrep'] in sym:
                        occ = mos.occs[iorb]
                        if abs(occ-2.) < 1.e-4:
                            occmap.append(iorb)
                        elif abs(occ) < 1.e-4:
                            virtmap.append(iorb)
                        else:
                            print(" Error: invalid occupation!", occ)
                            exit(5)

                nocc=len(occmap)
                nvirt=len(virtmap)

                print("\n Considering state: %s"%state['name'])
                print("  Number of occupied orbitals:", nocc)
                print("  Number of virtual orbitals:", nvirt)
                print("  Tensor space dimension:", nocc*nvirt)

                if os.path.exists('sing_%s'%state['irrep']):
                    readf = 'sing_%s'%state['irrep']
                    print('  Reading information of singlet calculation from file %s'%readf)
                elif os.path.exists('trip_%s'%state['irrep']):
                    readf = 'trip_%s'%state['irrep']
                    print('  Reading information of triplet calculation from file %s'%readf)
                elif self.ioptions['TDA'] and os.path.exists('ciss_%s'%state['irrep']):
                    readf = 'ciss_%s'%state['irrep']
                    print('  Reading information of TDA singlet calculation from file %s'%readf)
                elif self.ioptions['TDA'] and os.path.exists('cist_%s'%state['irrep']):
                    readf = 'cist_%s'%state['irrep']
                    print('  Reading information of TDA triplet calculation from file %s'%readf)
                else:
                    print('No file with information about the excited state (sing_a, trip_a, ...) found!')
                    exit(7)

                # the file is parsed once for every state
                # this is not efficient but also not critical
                # Note that this routine for non-TDA only reads X vector!
                curr_state = 0
                for line in open(readf):
                    if 'tensor space dimension' in line:
                        words = line.split()
                        space_dim = int(words[-1])
                        assert(nocc*nvirt == space_dim)
                    elif 'eigenvalue' in line:
                        words = line.split()
                        curr_state = float(words[0])
                        iiocc  = 0
                        iivirt = 0
                    elif curr_state == state['state_ind']:
                        words = [line[0+20*i:20+20*i] for i in range(4)]
                        for word in words:
                            state['tden'][occmap[iiocc],virtmap[iivirt]] = float(word.replace('D','E'))
                            iivirt+=1
                            if iivirt == nvirt:
                                iivirt = 0
                                iiocc += 1
                                if iiocc == nocc:
                                    curr_state = 0
                                    break
                    elif curr_state > state['state_ind']: break
 
        else: #UKS/UHF
            alpha_channel=(self.ioptions['spin'] >= 0)
            nocc={}
            nvirt={}
            for state in state_list:
                state['name'] = '%i%s'%(state['state_ind'],state['irrep'])
                state['tden'] = self.init_den(mos, rect=True)
                occmap  = []
                virtmap = []
                for iorb, sym in enumerate(mos.syms):
                    if state['irrep'] in sym:
                        occ = mos.occs[iorb]
                        if abs(occ-1.) < 1.e-4:
                            occmap.append(iorb)
                        elif abs(occ) < 1.e-4:
                            virtmap.append(iorb)
                        else:
                            print(" Error: invalid occupation!", occ)
                            exit(5)
                nocc=len(occmap)
                nvirt=len(virtmap)

                print("\n Considering state: %s"%state['name'])
                if alpha_channel:
                    print("  Number of occupied alpha orbitals:", nocc)
                    print("  Number of virtual alpha orbitals:", nvirt)
                    print("  Tensor space dimension (alpha):", nocc*nvirt)
                    assert(act_occ[0]==nocc)
                    assert(act_virt[0]==nvirt)
                else:
                    print("  Number of occupied beta orbitals:", nocc)
                    print("  Number of virtual beta orbitals:", nvirt)
                    print("  Tensor space dimension (beta):", nocc*nvirt)
                    assert(act_occ[1]==nocc)
                    assert(act_virt[1]==nvirt)

                if os.path.exists('unrs_%s'%state['irrep']):
                    readf = 'unrs_%s'%state['irrep']
                    print('  Reading information of unrestricted calculation from file %s'%readf)
                elif self.ioptions['TDA'] and os.path.exists('ucis_%s'%state['irrep']):
                    readf = 'ucis_%s'%state['irrep']
                    print('  Reading information of TDA unrestricted calculation from file %s'%readf)
                else:
                    print('No file with information about the excited state (sing_a, trip_a, ...) found!')
                    exit(7)

                curr_state = 0
                itdm = 0  # global index across alpha+beta blocks
                n_a = act_occ[0] * act_virt[0] #alpha space dimension (from turbomole directly)
                n_b = act_occ[1] * act_virt[1] #beta space dimension (from turbomole directly)
                n_total = n_a + n_b
                for line in open(readf):
                    if 'tensor space dimension' in line:
                        words = line.split()
                        space_dim = int(words[-1])
                        assert(n_total == space_dim)
                    elif 'eigenvalue' in line:
                        words = line.split()
                        curr_state = float(words[0])
                        iiocc = 0
                        iivirt = 0
                        itdm = 0  # reset global counter for new eigenstate
                    elif curr_state == state['state_ind']:
                        words = [line[20*i:20*(i+1)] for i in range(4)]
                        for word in words:
                            val = float(word.replace('D', 'E'))
                            #If TDA, unrs_a etc. store only X vectors
                            if self.ioptions['TDA']==True:
                                # only process the block we actually want
                                if alpha_channel and itdm >= n_a:
                                    itdm += 1
                                    continue   #we skip the beta block since we are running alpha channel
                                if not alpha_channel and itdm < n_a:
                                    itdm += 1
                                    continue   #we skip the alpha block since we are running beta channel
                            else:
                                if alpha_channel and itdm >= 2*n_a: #double size cause both X and Y are printed out in unrs_
                                    itdm += 1
                                    continue   #we skip the beta block since we are running alpha channel
                                if not alpha_channel and itdm < 2*n_a:
                                    itdm += 1
                                    continue   #we skip the alpha block since we are running beta channel
                            #This only reads X vector anyway, as done for RKS.
                            state['tden'][occmap[iiocc], virtmap[iivirt]] = val
                            iivirt += 1
                            if iivirt == nvirt:
                                iivirt = 0
                                iiocc += 1
                                if iiocc == nocc:
                                    curr_state = 0
                                    break
                            itdm += 1
                    elif curr_state > state['state_ind']:
                        break

        return state_list

    def ret_conf_tddft(self, rfile):

        # #Checking first which TURBOMOLE version is used. Output of MOs info changed from V8.0 unfortunately!
        # header=open(rfile,'r').readlines()[0:20]
        # version = None
        # for line in header:
        #     m = re.search(r'TURBOMOLE rev\.\s*V(\d+)-', line)
        #     if m:
        #         version = int(m.group(1))
        #         break

        rlines = open(rfile, 'r').readlines()[100:]
        ret_list = []
        occ_orb = False # section of the file
        restr = True #Default is restricted (RKS/RHF) calculations. Check in  turbomole output later
        act_occ=[] #List for active occupied out of turbomole output
        act_virt=[] #List for virtual occupied out of turbomole output

        for nr,line in enumerate(rlines):
            #Checking for UKS/UHF
            if 'CI SINGLES UHF-CALCULATION (spin-conserved)' in line or 'RPA UHF-EXCITATION-CALCULATION (spin-conserved)' in line:
                restr=False
                print("Detected unrestricted calculation. unrs_a or ucis_a files will be searched for.")
            # Checking for TDA
            if '(TAMM-DANCOFF-APPROXIMATION)' in line:
                self.ioptions['TDA']=True
                print("Detected TDA approximation.")
            # Reading n. of occupied and virtual orbitals for UKS or RKS. 
            # if version > 7: #Older TURBOMOLE versions have different output for MOs summary
            #     if restr is False and 'orbitals in total:  alpha spin   beta spin' in line:
            #         act_occ.extend(map(int, rlines[nr+3].split()[-2:]))
            #         act_virt.extend(map(int, rlines[nr+4].split()[-2:]))
            #     elif restr is True and 'orbitals in total:' in line:
            #         act_occ.append(int(rlines[nr+3].split()[-1]))
            #         act_virt.append(int(rlines[nr+4].split()[-1]))
            # else:
            if 'number of non-frozen orbitals' in line:
                    tmp_orb=int(line.split()[-1]) #Total n. of active orbitals (either alpha or beta channel)
                    tmp_occ=int(rlines[nr+1].split()[-1]) # Total n. of acive occupied orbitals (either alpha or beta channel)
                    act_occ.append(tmp_occ)
                    act_virt.append(tmp_orb-tmp_occ)
            #Extracting excitations
            if 'excitation' in line and not 'vector' in line and not 'energies' in line:
                ret_list.append({})
                words = line.split()
                ret_list[-1]['state_ind'] = int(words[0])
                #ret_list[-1]['irrep'] = file_handler.line_to_words(line)[2]
                if restr:
                    ret_list[-1]['irrep'] = words[2]
                else:
                    ret_list[-1]['irrep'] = words[1]
                ret_list[-1]['tot_en'] = eval(rlines[nr+3][40:])
                ret_list[-1]['exc_en'] = eval(rlines[nr+7][40:])
                ret_list[-1]['osc_str'] = eval(rlines[nr+16][40:])
                ret_list[-1]['char'] = []
            elif 'occ. orbital' in line:
                occ_orb = True
                #print 'occ. orbital in line'
            elif occ_orb:
                words = line.split()
                #print words
                if len(words) > 0:
                    #print 'words > 0'
                    if restr:
                        ret_list[-1]['char'].append({})
                        ret_list[-1]['char'][-1]['occ'] = words[0]+words[1]
                        ret_list[-1]['char'][-1]['virt'] = words[3]+words[4]
                        ret_list[-1]['char'][-1]['weight'] = float(words[-1])/100.
                        # escf does not print the phase of the coefficient
                    else:
                        ret_list[-1]['char'].append({})
                        ret_list[-1]['char'][-1]['occ'] = words[0] + words[1]
                        ret_list[-1]['char'][-1]['occ_spin'] = words[2]
                        ret_list[-1]['char'][-1]['virt'] = words[4] + words[5]
                        ret_list[-1]['char'][-1]['virt_spin'] = words[6]
                        ret_list[-1]['char'][-1]['weight'] = float(words[-1]) / 100.
                else:
                    occ_orb = False

        return ret_list,restr,act_occ,act_virt

#---

class file_parser_libwfa(file_parser_base):
    def __init__(self, ioptions):
        self.ioptions = ioptions
        self.excD = False
        self.excT = False
        self.state_list_om = [] # state_list read from the file ctnum*.om

    def read(self):
        return self.rmatfile_one()

    def rmatfile_one(self):
        """
        Read all om_at data from one file.
        """
        state_list = []

        mom = os.path.exists("ctnum_mulliken.om")
        lom = os.path.exists("ctnum_lowdin.om")

        if mom:
            if lom:
                if self.ioptions['Om_formula'] == 1:
                    fname = "ctnum_mulliken.om"
                elif self.ioptions['Om_formula'] == 2:
                    fname = "ctnum_lowdin.om"
                else:
                    raise error_handler.ElseError(str(self.ioptions['Om_formula']), "Om_formula")
            else:
                fname = "ctnum_mulliken.om"
        else:
            if lom:
                fname = "ctnum_lowdin.om"
            else:
                print("WARNING: Using old mode of reading .om files")
                return self.read_old()

        print("Reading %s ..."%fname)

        rfile=open(fname,'r')

        while True:
            try:
                line = next(rfile)
            except StopIteration:
                print("Finished parsing %s."%fname)
                break
            if len(line) < 3: continue

            # Parse the header
            state_list.append({})
            state = state_list[-1]
            if '<-->' in line: # ccman2 output
                line2 = line.split('<-->')[1]
                words = line2.split()
                print(line2)
                print(words)
            elif 'rho_imag' in line: # ccman2 output, complex
                l2 = line.split()[-1]
                words = [l2.split('_')[0]]
            else:
                words = line.split()
            if words[0] == 'Excited' and words[1] == 'State':
                state['name'] = 'es_%i'%int(words[2])
                state['exc_en'] = 0.
                state['osc_str'] = -1.
            elif words[0] == 'Spin-flipped' and words[1] == 'State':
                state['name'] = 'sf_%i'%int(words[2])
                state['exc_en'] = 0.
                state['osc_str'] = -1.
            else:
                state['name'] = words[0].replace('/', '-')
                state['exc_en'] = float(words[1]) * units.energy['eV'] if len(words) >= 2 else  0.
                state['osc_str'] = float(words[2]) if len(words) >= 3 else -1.
            state['lname'] = line.strip()

            # Dimensions
            words = next(rfile).split()
            dima=int(words[1])
            dimb=int(words[2])

            # Parse the array
            outarr = numpy.zeros([dima,dimb], float)

            ia = ib = 0
            while ia < dima:
              for val in next(rfile).split():
                outarr[ib, ia] = float(val)
                ib += 1
                if ib == dimb:
                  ia += 1
                  ib = 0
            state['OmAt'] = outarr
            state['Om'] = outarr.sum()

        rfile.close()
        return state_list

# Everything below here can be deleted at some point

    def read_old(self):
        """
        Read the ctnum*.om file as created with libwfa.
        This assumes the new output format where only one file is present.
        """
        state_list = []

        basedir='.'
        #ist = 1
        for omfile in sorted(os.listdir(basedir)):
            suff = omfile.split('.')[-1]
            if suff!='om': continue

            (typ, exctmp, osc, num_at, num_at1, om_at) = self.rmatfile(os.path.join(basedir,omfile))
            if typ == None: continue

            state_list.append({})

            state_list[-1]['lname'] = typ
            # Here several try statements are used to check for different versions of the output
            #    Not very elegant but it should work ...
            try:
                multl, iirrep, ist = typ.split('_')
            except:
                state_list[-1]['name'] = typ
            else:
                try:
                    state_list[-1]['irrep'] = self.ioptions.get('irrep_labels')[int(iirrep)]
                except:
                    state_list[-1]['irrep'] = iirrep

                state_list[-1]['mult'] = {'singlet':'(1)',
                                          'triplet':'(3)',
                                          'any':'(-)'}[multl]

                state_list[-1]['state_ind'] = int(ist)

                state_list[-1]['name'] = '%i%s%s'%(state_list[-1]['state_ind'],
                                                    state_list[-1]['mult'],
                                                    state_list[-1]['irrep'])

            #state_list[-1]['char'] = []
            state_list[-1]['exc_en'] = exctmp * units.energy['eV']
            state_list[-1]['osc_str'] = osc

            state_list[-1]['Om']   = om_at.sum()
            state_list[-1]['OmAt'] = om_at

          #ist += 1

        return state_list

    def rmatfile(self,fname):
        """
        Read Omega (atoms) matrix as produced by libwfa
        """
        print("Reading: %s ..."%fname)

        try:
            rfile=open(fname,'r')
        except IOError:
            print("\n WARNING: could not open %s."%fname)
            print("Did the calculation converge?\n")
            return None, None, None, None, None, None

        line=rfile.readline()
        words=line.split()

        if '<-->' in line:
            # ccman2 output
            typ = words[2]
            excen = float(words[3]) if len(words) >= 4 else  0.
            osc   = float(words[4]) if len(words) >= 5 else -1.
        else:
            typ = words[0]
            excen = float(words[1]) if len(words) >= 2 else  0.
            osc   = float(words[2]) if len(words) >= 3 else -1.

        line=rfile.readline()
        words=line.split()
        dima=int(words[1])
        dimb=int(words[2])

        #print " Dimensions: %i x %i"%(dima, dimb)

        outarr = numpy.zeros([dima,dimb])

        ia = 0
        ib = 0
        while 1:
          line = rfile.readline()
          if not line:
            break
          vals = [float(word) for word in line.split()]

          for val in vals:
            outarr[ib, ia] = val
            ib += 1
            if ib == dimb:
              ia += 1
              ib = 0
              if ia == dima:
                print(" All values read in")

        return typ, excen, osc, dima, dimb, outarr

    def parse_line(self, state, line, rfile=None):
        if 'Exciton analysis of the difference density matrix' in line:
            self.excD = True
            self.excT = False
        elif 'Exciton analysis of the transition density matrix' in line:
            self.excD = False
            self.excT = True

        self.parse_keys(state, self.excD, self.excT, line, rfile)

    def parse_keys(self, state, exc_diff, exc_1TDM, line, rfile=None):
        self.parse_key(state, 'osc_str', line, 'Osc. strength:')
        self.parse_key(state, 'mu', line, ' Dipole moment [D]', rfile=rfile)
        self.parse_key(state, 'Tmu', line, 'Trans. dipole moment [D]', rfile=rfile)
        self.parse_key(state, 'r2x', line, '  <r^2>', -3, not_string='Total')
        self.parse_key(state, 'r2y', line, '  <r^2>', -2, not_string='Total')
        self.parse_key(state, 'r2z', line, '  <r^2>', -1, not_string='Total')
        self.parse_key(state, 'r2', line, 'Total <r^2>')
        self.parse_key(state, 'nu', line, 'Number of unpaired electrons:', 2)
        self.parse_key(state, 'nunl', line, 'Number of unpaired electrons')
        self.parse_key(state, 'p', line, 'Number of detached / attached electrons')
        self.parse_key(state, 'Om_', line, 'omega', 2)
        self.parse_key(state, 'LOC', line, 'LOC', 2, not_string='LOCa')
        self.parse_key(state, 'LOCa', line, 'LOCa', 2)
        self.parse_key(state, 'QT2', line, 'QT2')
        self.parse_key(state, 'QTa', line, 'QTa')
        self.parse_key(state, 'Phe', line, '<Phe>', 2)
        self.parse_key(state, 'PRNTO', line, 'PR_NTO')
        self.parse_key(state, 'PRD', line, 'PR_D', 6)
        self.parse_key(state, 'PRA', line, 'PR_A')
        self.parse_key(state, '2P', line, 'Two-photon absorption cross-section [a.u.]', rfile=rfile, rmatrix=True)
        self.parse_key(state, 'S_HE', line, 'Entanglement entropy')
        self.parse_key(state, 'Z_HE', line, 'Nr of entangled states')
        self.parse_key(state, 'sigR', line, 'RMS size of the density', rfile=rfile)

        if exc_diff:
            self.parse_key(state, 'sigD', line, 'Hole size', rfile=rfile)
            self.parse_key(state, 'sigA', line, 'Electron size', rfile=rfile)
            self.parse_key(state, 'dD-A', line, '|<r_e - r_h>|')

        if exc_1TDM:
            self.parse_key(state, 'muT',  line, 'Trans. dipole moment [D]:')
            self.parse_key(state, 'dexc', line, 'RMS electron-hole separation', rfile=rfile)
            self.parse_key(state, 'sigH', line, 'Hole size', rfile=rfile)
            self.parse_key(state, 'sigE', line, 'Electron size', rfile=rfile)
            self.parse_key(state, 'dH-E', line, '|<r_e - r_h>|')
            self.parse_key(state, 'COV', line, 'Covariance(r_h, r_e) [Ang^2]')
            self.parse_key(state, 'Corr', line, 'Correlation coefficient')
            self.parse_key(state, 'rHx', line, '<r_h>', -3)
            self.parse_key(state, 'rHy', line, '<r_h>', -2)
            self.parse_key(state, 'rHz', line, '<r_h>', -1)
            self.parse_key(state, 'rEx', line, '<r_e>', -3)
            self.parse_key(state, 'rEy', line, '<r_e>', -2)
            self.parse_key(state, 'rEz', line, '<r_e>', -1)

class file_parser_qcadc(file_parser_libwfa):
    """
    Parse information from qchem.out in addition to the .om file.
    """
    def read(self):
        state_list = []

        basedir='.'

        self.state_list_om = self.rmatfile_one()

        exc_diff = False
        exc_1TDM = False
        self.irrep_labels = None
        rf = open(self.ioptions.get('rfile'), 'r')
        for line in rf:
            words = line.split()
            if 'Irreducible representations in point group:' in line:
                self.irrep_labels = line.lstrip('Irreducible representations in point group:').split()

            elif ' Term symbol' in line:
                if self.irrep_labels == None:
                    print("\n   WARNING: irrep labels not found in %s"%self.ioptions.get('rfile'))
                    print("   Use 'adc_print 2' or enter them into %s"%self.ioptions.ifile)
                    print("   Using info from %s: irrep_labels="%self.ioptions.ifile, self.ioptions.get('irrep_labels'))
                    print()
                    self.irrep_labels = self.ioptions.get('irrep_labels')

                state_list.append({})
                state = state_list[-1]

                state['state_ind'] = int(words[2])
                state['mult']      = words[3]
                state['irrep']     = words[4]
                state['name']      = '%s%s%s'%(words[2], words[3], words[4])

                oname = self.om_file_name(state)
                for ostate in self.state_list_om:
                    if oname == ostate['lname'].split()[0]:
                        for prop in ['Om', 'OmAt', 'exc_en', 'osc_str']:
                            state[prop] = ostate[prop]

            elif 'MP(2) Summary' in line:
                state_list.append({})
                state_list[-1]['name'] = 'gr-st'
                state_list[-1]['exc_en'] = 0.
                state_list[-1]['state_ind'] = 0

            elif ' Excitation energy:' in line:
                exc_chk = float(words[2])

                if not 'exc_en' in state_list[-1]:
                    state_list[-1]['exc_en'] = exc_chk

                if abs(exc_chk - state_list[-1]['exc_en']) > exc_chk * 1.e-4:
                    print(exc_chk, state_list[-1]['exc_en'])
                    raise error_handler.MsgError("Excitation energies do not match\n   -> Try using rtype='libwfa' instead.")

            elif 'Exciton analysis of the difference density matrix' in line:
                exc_1TDM = False
                if len(state_list) > 0: exc_diff = True

            elif 'Exciton analysis of the transition density matrix' in line:
                exc_diff = False
                if len(state_list) > 0: exc_1TDM = True

            elif 'Transition Summary' in line:
                break

            if len(state_list) > 0:
                self.parse_keys(state_list[-1], exc_diff, exc_1TDM, line, rf)

        rf.close()

        # Post-processing
        # Pre-factor for 2P cross-section
        pre2P = numpy.pi**3 * units.constants['c0']**(-2) * units.tpa['GM']
        print("pre2P, GM: ", pre2P, units.tpa['GM'])
        # This is numerically the same as [JCP, 146, 174102]:
        #   numpy.pi**3 / u.c0 * u.cm**5 / 2.9979E10 * 10**50
        for state in state_list:
            # 2P cross-section pre-factor in GM
            if '2P' in state:
                state['GM'] = pre2P * state['2P']/30 * (state['exc_en'] / units.energy['eV'])
                state['2P'] *= .000001
            if 'Tmux' in state:
                state['Tmu'] = (state['Tmux']**2 + state['Tmuy']**2 + state['Tmuz']**2)**.5

        return state_list

    def om_file_name(self, state):
        """
        Construct the name of the state in the .om file
        """
        multo = {'(1)':'singlet',
                 '(3)':'triplet',
                 '(-)':'any'}[state['mult']]

        # subtract 1 for ground state irrep
        state_indo = state['state_ind']
        if state['mult'] == '(1)':
            if state['irrep'].lower() in ['a', 'ag', 'a1', "a'"]:
                state_indo -= 1

        return '%s_%s_%i'%(multo, state['irrep'], state_indo)

class file_parser_qctddft(file_parser_libwfa):
    def read(self, mos):
        """
        Read the X vector from standard output. Y is discarded.
        """
        state_list = []
        exc_diff = exc_1TDM = tdread = libwfa = False
        spin_flip = False
        istate = 0

        self.state_list_om = self.rmatfile_one()

        if self.ioptions.get('TDA'):
            ststr  = 'TDDFT/TDA Excitation Energies'
            ststr2 = 'CIS Excitation Energies'
        else:
            ststr  = 'TDDFT Excitation Energies'
            ststr2 = 'xyzabc' # dummy string

        print("Parsing %s for %s ..."%(self.ioptions.get('rfile'), ststr))

        nsing = 0
        ntrip = 0
        rfileh = open(self.ioptions.get('rfile'),'r')
        while True: # loop over all lines
            try:
                line = next(rfileh)
            except StopIteration:
              print("Finished parsing file %s"%self.ioptions.get('rfile'))
              break

            if ststr in line or ststr2 in line:
                tdread = True
                line = next(rfileh)
            elif 'SF-DFT Excitation Energies' in line:
                # SF seems to be only possible with TDA
                # -> no need to check for TDA/RPA
                tdread = True
                line = next(rfileh)
                line = next(rfileh)
                line = next(rfileh)
                spin_flip = True
            elif 'TDDFT calculation will be performed' in line:
                tdread = False
            elif '-----' in line:
                tdread = False
            elif 'Excited State Analysis' in line:
                tdread = False
                libwfa = True
                if len(state_list) == 0:
                    errstr = "No excitation energy parsed!"
                    errstr+= "\n   Please, set 'TDA=True' if this was a TDA calculation."
                    raise error_handler.MsgError(errstr)
            elif 'SA-NTO Decomposition' in line:
                libwfa = False
            elif 'Welcome to Q-Chem' in line:
                if len(state_list) > 0:
                    print("\n WARNING: found second Q-Chem job!\n   Deleting everything parsed so far.\n")

                state_list = []
                exc_diff = exc_1TDM = tdread = libwfa = False
                istate = 0

            if tdread:
                words = line.replace(':','').split()
                if 'Excited state' in line:
                    state_list.append({})
                    state = state_list[-1]

                    state['state_num'] = int(words[2])
                    state['exc_en'] = float(words[-1])

                    line = next(rfileh)
                    line = next(rfileh)
                    if 'J1' in line:
                        line = self.read_ecomp(line, rfileh, state)
                    words = line.split()

                    if words[0] == 'Multiplicity:':
                        state['mult'] = words[1]
                    else:
                        state['mult'] = 'X'

                    if state['mult'] == 'Singlet':
                        nsing += 1
                        state['name'] = "S_%i"%nsing
                        state['lname'] = "Singlet %i"%nsing
                    elif state['mult'] == 'Triplet':
                        ntrip += 1
                        state['name'] = "T_%i"%ntrip
                        state['lname'] = "Triplet %i"%ntrip
                    elif spin_flip:
                        state['name'] = 'sf_%i'%(state['state_num'])
                        state['lname'] = 'Spin-flipped State %i'%(state['state_num'])
                    else:
                        state['name'] = 'es_%i'%(state['state_num'])
                        state['lname'] = 'Excited State %i'%(state['state_num'])

                    if self.ioptions['read_libwfa']:
                        om_at = None
                        for ostate in self.state_list_om:
                            if ostate['lname'] == state['lname']:
                                om_at = ostate['OmAt']

                        if om_at is None:
                            if os.path.exists('%s_ctnum_atomic.om'%state['name']):
                                (typ, exctmp, osc, num_at, num_at1, om_at) = self.rmatfile('%s_ctnum_atomic.om'%state['name'])
                            elif os.path.exists('%s_ctnum_mulliken.om'%state['name']):
                                (typ, exctmp, osc, num_at, num_at1, om_at) = self.rmatfile('%s_ctnum_mulliken.om'%state['name'])
                            elif os.path.exists('%s_ctnum_lowdin.om'%state['name']):
                                (typ, exctmp, osc, num_at, num_at1, om_at) = self.rmatfile('%s_ctnum_lowdin.om'%state['name'])

                        if not om_at is None:
                            state_list[-1]['Om']   = om_at.sum()
                            state_list[-1]['OmAt'] = om_at
                    else:
                        state_list[-1]['tden'] = self.init_den(mos, rect=True)

                elif 'Strength' in line:
                    state_list[-1]['osc_str'] = float(words[-1])

                elif 'amplitude' in line and not self.ioptions['read_libwfa']:
                    # ignore the Y vector.
                    #    Otherwise the Y would go into the virt-occ block!
                    if 'Y:' in line: continue

                    awords = self.delete_chars(line, ['X:', 'Y:', 'D', 'V', '(', ')', '-->'])

                    iocc = int(awords[0]) - 1
                    ivirt = int(awords[1]) + mos.ret_ihomo()

                    coeff =  float(awords[4])

                    state_list[-1]['tden'][iocc, ivirt] += coeff

                elif 'J1' in line:
                    line = self.read_ecomp(line, rfileh, state_list[-1])

            if libwfa:
                if 'Excited state' in line:
                    words = line.replace(':', '').split()
                    istate = int(words[2]) - 1

                # Disentangle the order of singlets and triplets
                elif '  Singlet' in line or '  Triplet' in line:
                    words = line.split()
                    currname = words[0][0] + "_" + words[1]
                    for istate, state in enumerate(state_list):
                        if state['name'] == currname:
                            break
                    else:
                        raise error_handler.MsgError("Did not find state %s"%currname)

                elif 'Exciton analysis of the difference density matrix' in line:
                    exc_1TDM = False
                    exc_diff = True

                elif 'Exciton analysis of the transition density matrix' in line:
                    exc_diff = False
                    exc_1TDM = True

                self.parse_keys(state_list[istate], exc_diff, exc_1TDM, line)

        rfileh.close()
        return state_list

    def read_ecomp(self, line, rfileh, state):
        """
        Read energy components if available.
        """
        words = line.split()
        state['F1']  = float(words[1]) + float(words[3]) + float(words[5]) + float(words[7])
        state['F1K'] = float(words[1]) + float(words[3]) + float(words[7])
        state['K1']  = float(words[5])
        line = next(rfileh)
        words = line.split()
        state['J2'] = float(words[1])
        state['K2'] = float(words[3])
        state['K12']= state['K1'] + state['K2']
        state['XC2'] = float(words[5])
        line = next(rfileh)
        if 'J2' in line:
            words = line.split()
            state['J2'] = float(words[2])
            line = next(rfileh)

        return line

#---

class file_parser_col(file_parser_base):
    def read_iwfmt(self, dens, filen, fac = 1.):
        """
        Read output from iwfmt.x for a 1-particle density file.
        """
        eref = eexc = 0.

        in_data=False
        #for line in open(filen,'r').readlines():
        rfile = open(filen, 'r')
        while True:
            try:
                line=next(rfile)
            except StopIteration:
                break
            if ' 0.000000000000E+00' in line:
                if len(line.split()) == 1:
                    words = last_line.split()
                    (num,lab1,ibvtyp,itypea,itypeb,ifmt,last,nipv) = [int(word) for word in words]
                    if itypea==0 and itypeb==7:
                        print(' Reading symmetric density ...')
                        in_data=True
                        sym = 1
                    elif itypea==2 and itypeb==9:
                        print(' Reading antisymmetric density ...')
                        in_data=True
                        sym = -1
                    else:
                        print('Warning: this section does not contain a density')
                        print('itypea: %i, itypeb: %i'%(itypea,itypeb))
                        print(last_line)
                        in_data=False
            elif '-1025   -1025      -1' in line:
                line=next(rfile)
                (eexc, eref, ecore) = [float(word) for word in line.split()]
            elif '    ' in line:
                pass
            elif in_data:
                words = line.split()

                val = float(words[0])
                i = int(words[1])-1
                j = int(words[2])-1

                if 0.2 < abs(val) < 1.9999: print("(i,j)=(%2i,%2i), val=%6.3f"%(i,j,val))
                dens[j,i] += fac * val

                if i != j: dens[i,j] += fac * sym * val

            last_line = line

        return eref, eexc

    def read_sifs(self):
        """
        Directly read a SIFS file.
        This would require a small Fortran program with the required routines.
            -> Should be compiled with Columbus.
        Interface with f2py
        """
        raise error_handler.NIError()

    def read_trncils(self, state, mos, filen):
        """
        Read output from transci.x for a 1-particle density file.
        """
        # change the format of the MO labels
        mos.syms2 = [sym.lower().replace('_','') for sym in mos.syms]
        if not self.ioptions['ncore'] == {}:
            for isym, sym2 in enumerate(mos.syms2):
                imo, irr = self.sym_split(sym2)
                new = "%i%s"%(imo-self.ioptions['ncore'][irr], irr)
                mos.syms2[isym] = new

        tmp = filen.replace('LISTINGS/','').replace('state','').replace('drt','').replace('trncils.FROM','').split('TO')
        tmp2 = tmp[0].split('.')
        ir_st_ref = (int(tmp2[0]), int(tmp2[1]))
        tmp2 = tmp[1].split('.')
        ir_st_exc = (int(tmp2[0]), int(tmp2[1]))

        state['state_ind'] = int(ir_st_exc[1])
        #state['name'] = '%s.%i-%s.%i'%(irrep_ref, int(ir_st_ref[1]), state['irrep'], state['state_ind'])
        state['name'] = 'DRT%i.St%02i'%ir_st_exc

        if ir_st_ref==ir_st_exc:
            state['sden'] = self.init_den(mos)
            dfac = 1.
            mat = state['sden']
        else:
            state['tden'] = self.init_den(mos)
            dfac = 1. / numpy.sqrt(2)
            mat = state['tden']

        rfile = open(filen,'r')
        while True:
            try:
                line=next(rfile)
            except StopIteration:
                break
            if ' symm-trans density matrix:' in line:
                print('file %s, reading:'%filen)
                print(line.strip('\n'))
                line=next(rfile)
                self.read_block_mat(mat, mos, rfile, sym=1, dfac=dfac)
            elif 'asymm-trans density matrix:' in line:
                print('reading ...')
                print(line[:-1])
                line=next(rfile)
                self.read_block_mat(mat, mos, rfile, sym=-1, dfac=dfac)
            elif 'Transition energy:  ' in line:
                words=line.split() # use the line with "a.u." because it is given with more digits
                state['exc_en']=float(words[2]) * units.energy['eV']
                next(rfile)
            elif 'Oscillator strength : ' in line:
                words = line.split()
                state['osc_str']=float(words[-1])
        rfile.close()

        try:
            self.read_ciudgsm(state, ir_st_ref, ir_st_exc)
        except:
            print("   WARING: could not read extensivity corrected energies")

    def read_ciudgsm(self, state, ir_st_ref, ir_st_exc):
        """
        Read extensivity corrected energies from ciudgsm.
        """
        filen_ref = 'LISTINGS/ciudgsm.drt%i.sp'%ir_st_ref[0]
        filen_exc = 'LISTINGS/ciudgsm.drt%i.sp'%ir_st_exc[0]

        eci_ref = []
        print('Reading %s ...'%filen_ref)
        for line in open(filen_ref, 'r'):
            if 'eci' in line:
                words = line.split()
                eci_ref.append(float(words[2]))

        if ir_st_exc[0] == ir_st_ref[0]:
            eci_exc = eci_ref
        else:
            print('Reading %s ...'%filen_exc)
            eci_exc = []
            for line in open(filen_exc, 'r'):
                if 'eci' in line:
                    words = line.split()
                    eci_exc.append(float(words[2]))

        iref, iexc = ir_st_ref[1] - 1, ir_st_exc[1] - 1

        tmp_en = (eci_exc[5*iexc] - eci_ref[5*iref]) * units.energy['eV']
        assert(abs(state['exc_en'] - tmp_en) < 1.e-6)
        state['+DV1'] = (eci_exc[5*iexc+1] - eci_ref[5*iref+1]) * units.energy['eV']
        state['+DV2'] = (eci_exc[5*iexc+2] - eci_ref[5*iref+2]) * units.energy['eV']
        state['+DV3'] = (eci_exc[5*iexc+3] - eci_ref[5*iref+3]) * units.energy['eV']
        state['+P']   = (eci_exc[5*iexc+4] - eci_ref[5*iref+4]) * units.energy['eV']

    def read_block_mat(self, mat, mos, rfile, sym, dfac=1.):
        """
        Parse the block matrix output in the listing file.
        This is not the cleanest routine but it seems to work ...
        """
        while(1):
            words=next(rfile).replace('MO','').split()

            if 'density' in words: # skip two lines in case the next symmetry block is coming
                line=next(rfile)
                words=next(rfile).replace('MO','').split()

            try:
                head_inds = [mos.syms2.index(sym2) for sym2 in words]
            except:
                print(" ERROR reading:", words)
                print("syms2: ", mos.syms2)
                print()
                raise

            words=next(rfile).replace('MO','').split()
            while(len(words)>0): # loop over a block with constant header labels
                if 'integral' in words: break
                try:
                    left_ind = mos.syms2.index(words[0])
                except:
                    print('\n ERROR parsing: ')
                    print(words)
                    raise

                for i, word in enumerate(words[1:]):
                    val = float(word)
                    if val*val > 0.01: print("(%2i->%2i)-(%s->%s), val=% 8.4f"%\
                       (head_inds[i],left_ind,mos.syms2[head_inds[i]],mos.syms2[left_ind],val))

                    mat[head_inds[i],left_ind] += dfac * val
                    if left_ind != head_inds[i]:
                        mat[left_ind,head_inds[i]] += dfac * sym*val

                words=next(rfile).replace('MO','').split()

            if 'integral' in words: break

class file_parser_col_mrci(file_parser_col):
    def read(self, mos):
        if self.ioptions['s_or_t'] == 's':
            raise error_handler.MsgError('analyze_sden.py not implemented for col_mrci! Use "nos" instead.')
            # it is actually possible to have bra=ket transitions in transmomin ...

        state_list = []

        for lfile in sorted(os.listdir('LISTINGS')):
            if not 'trncils' in lfile: continue

            tmp  = lfile.replace('state','').replace('drt','').replace('trncils.FROM','').split('TO')
            tmp2 = tmp[0].split('.')
            ir_st_ref = (int(tmp2[0]), int(tmp2[1]))
            if ir_st_ref == self.ioptions['ref_state']:
                print("Reading %s ..."%lfile)
                state = {}
                self.read_trncils(state, mos, 'LISTINGS/%s'%lfile)
                if 'tden' in state:
                    state_list.append(state)
            else:
                print("Skipping %s (does not match ref_state)"%lfile)
                print(ir_st_ref)
                print(self.ioptions['ref_state'])

        return state_list

    def read_fcd(self, mos):
        state_list = [{}, {}, {}]

        ist1 = self.ioptions['state_pair'][0]
        ist2 = self.ioptions['state_pair'][1]

        fn1  = "LISTINGS/trncils.FROMdrt1.state%iTOdrt1.state%i"%(ist1, ist1)
        fn2  = "LISTINGS/trncils.FROMdrt1.state%iTOdrt1.state%i"%(ist2, ist2)
        fn12 = "LISTINGS/trncils.FROMdrt1.state%iTOdrt1.state%i"%(ist1, ist2)
        fn21 = "LISTINGS/trncils.FROMdrt1.state%iTOdrt1.state%i"%(ist2, ist1)

        self.read_trncils(state_list[0], mos, fn1)
        self.read_trncils(state_list[1], mos, fn2)
        if os.path.exists(fn12):
            self.read_trncils(state_list[2], mos, fn12)
        elif os.path.exists(fn21):
            self.read_trncils(state_list[2], mos, fn21)
        else:
            print("Not using transition charges.")
            state_list[2] == None

        return state_list

class file_parser_col_mcscf(file_parser_col):
    def read(self, mos):
        state_list = []

        for lfile in sorted(os.listdir('WORK')):
            # Find the suitable files. This could also be done with regexps ...
            if not '.iwfmt' in lfile: continue
            if not 'mcsd1fl' in lfile: continue

            if self.ioptions['s_or_t'] == 't':
                if not '-' in lfile: continue

                print("Reading %s ..."%lfile)
                state_list.append({})
                self.read_mc_tden(state_list[-1], mos, lfile)

            elif self.ioptions['s_or_t'] == 's':
                if '-' in lfile: continue

                print("Reading %s ..."%lfile)
                state_list.append({})
                self.read_mc_sden(state_list[-1], mos, lfile)

        if len(state_list) == 0:
            raise error_handler.MsgError('No density file found! Did you run write_den.bash?')

        return state_list

    def read_fcd(self, mos):
        state_list = [{}, {}, {}]

        ist1 = self.ioptions['state_pair'][0]
        ist2 = self.ioptions['state_pair'][1]

        fn1  = "mcsd1fl.drt1.st%2.2i.iwfmt"%ist1
        fn2  = "mcsd1fl.drt1.st%2.2i.iwfmt"%ist2
        fn12 = "mcsd1fl.drt1.st%2.2i-st%2.2i.iwfmt"%(ist1, ist2)
        fn21 = "mcsd1fl.drt1.st%2.2i-st%2.2i.iwfmt"%(ist2, ist1)

        self.read_mc_sden(state_list[0], mos, fn1)
        self.read_mc_sden(state_list[1], mos, fn2)
        try:
            self.read_mc_tden(state_list[2], mos, fn12)
        except IOError:
            print("File %s not found, trying %s..."%(fn12, fn21))
            try:
                self.read_mc_tden(state_list[2], mos, fn21)
            except IOError:
                print("Not using transition charges.")
                state_list[2] == None

        return state_list

    def read_mc_tden(self, state, mos, filen):
        tmp = filen.replace('mcsd1fl.drt','').replace('st','').replace('.iwfmt','').split('.')
        state['irrep'] = self.ioptions.get('irrep_labels')[int(tmp[0]) - 1]
        state['name'] = '%s.%s'%(state['irrep'], tmp[1])
        state['tden'] = self.init_den(mos)

        # read symmetric part
        (eref, eexc) = self.read_iwfmt(state['tden'], 'WORK/'+filen, fac=1/numpy.sqrt(2.))
        state['exc_en'] = (eexc - eref) * units.energy['eV']

        # read antisymmetric part
        self.read_iwfmt(state['tden'], 'WORK/'+filen.replace('mcsd1fl', 'mcad1fl'), fac=1/numpy.sqrt(2.))

    def read_mc_sden(self, state, mos, filen):
        tmp = filen.replace('mcsd1fl.drt','').replace('st','').replace('.iwfmt','').split('.')
        state['state_ind'] = int(tmp[1])
        state['name'] = 'DRT%i.St%02i'%(int(tmp[0]), int(tmp[1]))
        state['sden'] = self.init_den(mos)

        (eref, eexc) = self.read_iwfmt(state['sden'], 'WORK/'+filen)
        state['exc_en'] = eexc

class file_parser_terachem(file_parser_base):
    """
    Read terachem TDDFT job.
    """
    def read(self, mos):
        state_list = []

        rfileh = open(self.ioptions['rfile'], 'r')
        while True: # loop over all lines
            try:
                line = next(rfileh)
            except StopIteration:
              print("Reached end of file %s"%self.ioptions.get('rfile'))
              break

            if 'Largest CI coefficients' in line:
                state_list.append({})
                state = state_list[-1]
                state['tden'] = self.init_den(mos, rect=True)

                while True:
                    line = next(rfileh)
                    words = line.split()
                    if len(words) == 0: break

                    iocc  = int(words[0]) - 1
                    avirt = int(words[2]) - 1
                    if ('X' in line) and ('Y' in line): # rpa case
                        coeff = float(line.split(':')[2].split('Y')[0])
                    else:
                        coeff = float(words[-1])

                    state['tden'][iocc, avirt] = coeff

            elif 'Final Excited State Results' in line:
                line = next(rfileh)
                line = next(rfileh)
                line = next(rfileh)
                for state in state_list:
                    line = next(rfileh)
                    words = line.split()

                    state['state_ind'] = int(words[0])
                    state['name'] = 'A' + words[0]
                    state['exc_en'] = float(words[2])
                    state['osc_str'] = float(words[3])

        rfileh.close()

        return state_list

class file_parser_tddftb(file_parser_base):
    """
    Read DFTB+ TDDFTB job.
    Author: Ljiljana Stojanovic
    """
    def read(self, mos):
        state_list = []  # check if state_list is already allocated
        occ = []
        virt = []
        words = []
        Coeff = False
        # read the order of iocc, avirt pairs from the SPX.DAT file. It corresponds to the order in which X+Y coeffs are written in XplusY.DAT file
        spxfileh = open('SPX.DAT', 'r') # spxfile is SPX.DAT
        for line in spxfileh:
            words = line.split()
            if len(words) == 6:
               occ.append(int(words[3]))
               virt.append(int(words[5]))
        spxfileh.close()

        words = []
        rfileh = open('XplusY.DAT', 'r') # read X+Y coeffs from XplusY.DAT
        k = 0
        for line in rfileh:
            words = line.split()
            if 'S' in line:
                 k += 1
                 Coeff = True
                 xply = []
                 state_list.append({})
                 state = state_list[-1]
                 state['tden'] = self.init_den(mos, rect=True)
            elif Coeff:
               try:
                 if words[1] != 'S':
                     i = -1
                     for x in range (0,len(words)):
                         i += 1
                         xply.append(float(words[i]))
                 else: continue
                 j = -1
                 for x in range (0,len(xply)-1):
                     j += 1
                     iocc = occ[j]-1
                     avirt = virt[j]-1
                     state['tden'][iocc, avirt] = xply[j]
               except:
                 if len(words) == 3 and words[1]=='S': continue
               state['tden']=state['tden']/numpy.linalg.norm(state['tden'])

        rfileh.close()

        excfile = open(self.ioptions['rfile'], 'r') # rfile is EXC.DAT
        line = next(excfile)
        line = next(excfile)
        line = next(excfile)
        line = next(excfile)
        line = next(excfile)
        state_ind = 0
        for state in state_list:
            line = next(excfile)
            words = line.split()
            if len(words) == 8:
               state_ind += 1
               state['state_ind'] = state_ind
               state['name'] = str(state_ind) + 'A'
               state['exc_en'] = float(words[0])
               state['osc_str'] = float(words[1])

        excfile.close()

        return state_list

class file_parser_dftmrci(file_parser_base):
    """
    Read DFT/MRCI job.
    """
    def read(self, mos):
        state_list = []
        self.orb_map = numpy.zeros([mos.ret_num_mo(), mos.ret_num_mo()], int)

        self.read_stdout(mos)
        if self.ioptions.get('rfile2', False) == 'mrci.tcidens':
          self.read_tcidens(mos, state_list)
        else:
          self.read_cidens(mos, state_list)

        return state_list

    def read_stdout(self, mos):
        """
        Parse the standard output to get the mapping of MO indices with symmetry.
        """
        with open(self.ioptions['rfile'], 'r') as f:
            while True:
                try:
                    line = next(f)
                except StopIteration:
                    raise error_handler.MsgError('Did not find orbital map')

                if 'typ nr (CI/TM)' in line:
                    for iorb in range(mos.ret_num_mo()):
                        line = next(f)

                        # Frozen orbitals are kept as zero entries
                        if line.strip() == '.':
                            break
                        if line[:6] == 'frozen':
                            continue

                        words = line[7:].split()
                        ien = int(words[0])
                        isy = int(words[1])
                        self.orb_map[ien-1, isy-1] = 1
                    break

    def read_cidens(self, mos, state_list):
        lvprt = self.ioptions['lvprt']

        df = open('mrci.cidens', 'rb')
        print('\nReading file %s ...'%df.name)
        dummy, nroot = struct.unpack('2i', df.read(8))
        print('nroot:', nroot)

        en0 = None
        for iroot in range(2*nroot-1):
            df.read(8)
            lab1 = struct.unpack('6s', df.read(6))[0].strip().decode('UTF-8')
            lab2 = struct.unpack('6s', df.read(6))[0].strip().decode('UTF-8')
            df.read(8) # dummy
            en = struct.unpack('d', df.read(8))[0]

            df.read(8)
            nmo, dummy  = struct.unpack('2i', df.read(8))
            dummy, sden = struct.unpack('2i', df.read(8))
            dummy, reclen = struct.unpack('2i', df. read(8))

            if lvprt >= 2:
                print("\nParsing <%s|E_pq|%s>.  nmo: %i, reclen: %i"%(lab1, lab2, nmo, reclen))

            if nmo != mos.ret_num_mo():
                raise error_handler.MsgError('Inconsistent number of MOs')
            nentry = nmo * nmo
            if nentry*8 != reclen:
                raise error_handler.MsgError('Inconsistent record length')

            coeff = struct.unpack(nentry*'d', df.read(nentry*8))

            if self.ioptions['s_or_t'] == 's':
                if sden != 0:
                    print(" Analysing state %s"%lab1)
                    state_list.append({})
                    state = state_list[-1]
                    state['name'] = lab1.replace("'", "1").replace('"', "11")
                    if en0 is None:
                        state['exc_en'] = 0.
                        en0 = en
                    else:
                        state['exc_en'] = (en - en0) * units.energy['eV']

                    den_tmp = numpy.reshape(coeff, [nmo, nmo])
                    state['sden'] = numpy.dot(self.orb_map,\
                        numpy.dot(den_tmp, self.orb_map.T))
                    if lvprt >= 2:
                        self.dens_stat(state['sden'], lvprt)
            else:
                if sden == 0:
                    print(" Analysing transition %s -> %s"%(lab1, lab2))
                    state_list.append({})
                    state = state_list[-1]
                    state['name'] = lab2.replace("'", "1").replace('"', "11")
                    state['exc_en'] = en * units.energy['eV']

                    den_tmp = .5**.5 * numpy.reshape(coeff, [nmo, nmo])
                    state['tden'] = numpy.dot(self.orb_map,\
                        numpy.dot(den_tmp, self.orb_map.T))
                    if lvprt >= 2:
                        self.dens_stat(state['tden'], lvprt)

        df.close()

    def read_tcidens(self, mos, state_list):
        lvprt = self.ioptions['lvprt']

        df = open('mrci.tcidens', 'rb')
        print('\nReading file %s ...'%df.name)
        dummy, nmat = struct.unpack('2i', df.read(8))
        print('nmat:', nmat)

        en0 = None
        for imat in range(nmat):
            df.read(8)
            lab1 = struct.unpack('6s', df.read(6))[0].strip().decode('UTF-8')
            lab2 = struct.unpack('6s', df.read(6))[0].strip().decode('UTF-8')
            df.read(8) # dummy
            en = struct.unpack('d', df.read(8))[0]

            df.read(8)
            nmo, dummy  = struct.unpack('2i', df.read(8))
            dummy, sden = struct.unpack('2i', df.read(8))
            dummy, reclen = struct.unpack('2i', df.read(8))

            if lvprt >= 2:
                print("\nParsing <%s|E_pq|%s>.  nmo: %i, reclen: %i"%(lab1, lab2, nmo, reclen))

            if nmo != mos.ret_num_mo():
                raise error_handler.MsgError('Inconsistent number of MOs')
            nentry = nmo * nmo
            if nentry*8 != reclen:
                raise error_handler.MsgError('Inconsistent record length')

            coeff = struct.unpack(nentry*'d', df.read(nentry*8))

            print(" Analysing transition %s -> %s"%(lab1, lab2))
            state_list.append({})
            state = state_list[-1]
            state['name'] = lab2
            state['exc_en'] = en * units.energy['eV']

            den_tmp = .5**.5 * numpy.reshape(coeff, [nmo, nmo])
            state['tden'] = numpy.dot(self.orb_map,\
                numpy.dot(den_tmp, self.orb_map.T))
            if lvprt >= 2:
                self.dens_stat(state['tden'], lvprt)

        df.close()

class adf_kffile:
    """
    Wrapper for ADF.
    """
    def __init__(self, fname):
        try:
            from scm.plams import KFFile
        except ImportError:
            from kf import kffile as KFFile

        self.kffile = KFFile(fname)

    def read(self, sec, key):
        tmp_arr = numpy.array(self.kffile.read(sec, key))

        if len(tmp_arr.shape) == 0:
            return numpy.array([tmp_arr])
        else:
            return tmp_arr

class file_parser_adf(file_parser_base):
    """
    Read ADF TDDFT using the TAPE21 file.
    """
    def read(self, mos):
        rfile = adf_kffile(self.ioptions['rfile'])
        try:
            nmo = int(rfile.read('A','nmo_A'))
            nelec = int(rfile.read('General','electrons'))
        except:
            print(("\n  ERROR: reading TAPE21 file (%s)!\n"%self.ioptions['rfile']))
            raise
        assert nelec%2==0, "Odd number of electrons not supported"

        group = rfile.read('Geometry','grouplabel')[0].strip()
        if not group in ['NOSYM', 'C1', 'c1']:
            print(("grouplabel: " + group))
            raise error_handler.MsgError('No support for symmetry')

        nocc = nelec // 2
        assert nocc == mos.ret_ihomo() + 1
        nvirt = nmo - nocc
        assert nvirt == mos.ret_num_mo() - nocc

        try:
            self.nsing = len(rfile.read('All excitations','All Sing-Sing excitations'))
            excss = rfile.read('Excitations SS A', 'excenergies') * units.energy['eV']
            oscss = rfile.read('Excitations SS A', 'oscillator strengths')
        except:
            self.nsing = 0
        try:
            self.ntrip = len(rfile.read('All excitations','All Sing-Trip excitations'))
            excst = rfile.read('Excitations ST A', 'excenergies') * units.energy['eV']
        except:
            self.ntrip = 0

        state_list = [{} for istate in range(self.nsing+self.ntrip)]
        istate = 0
        for ising in range(self.nsing):
            state = state_list[istate]

            state['state_ind'] = ising + 1
            state['exc_en'] = excss[ising]
            state['osc_str'] = oscss[ising]
            state['mult'] = 1
            state['name'] = '%i(%i)A'%(state['state_ind'], state['mult'])

            state['tden'] = self.init_den(mos, rect=True)
            eigen = rfile.read('Excitations SS A','eigenvector %s'%(istate+1))
            state['tden'][:,nocc:nmo] = eigen.reshape(nocc, nvirt)
            istate += 1

            # print-out
            print((state['name']))
            tden = state['tden']
            for i in range(len(tden)):
                for j in range(len(tden[0])):
                    val = tden[i, j]
                    if val*val > 0.1:
                        print(("(%i -> %i) % .4f"%(i+1,j+1,val)))

        for itrip in range(self.ntrip):
            state = state_list[istate]

            state['state_ind'] = itrip + 1
            state['exc_en'] = excst[itrip]
            state['mult'] = 3
            state['name'] = '%i(%i)A'%(state['state_ind'], state['mult'])

            state['tden'] = self.init_den(mos, rect=True)
            eigen = rfile.read('Excitations ST A','eigenvector %s'%(istate+1 - self.nsing))
            state['tden'][:,nocc:nmo] = eigen.reshape(nocc, nvirt)
            istate += 1

            # print-out
            print((state['name']))
            tden = state['tden']
            for i in range(len(tden)):
                for j in range(len(tden[0])):
                    val = tden[i, j]
                    if val*val > 0.1:
                        print(("(%i -> %i) % .4f"%(i+1,j+1,val)))

        #rfile.close()

        return state_list

class file_parser_nos(file_parser_base):
    """
    Interpret the MO-file as a diagonal density.
    """
    def read(self, mos):
        state_list = []

        for no_file in self.ioptions['ana_files']:
            state_list.append({})
            state = state_list[-1]
            if not self.ioptions['unrestricted']:
                state['sden'] = self.read_no_file(state, mos, no_file)
            else:
                # TODO: One could compute the spin-density here
                state['sden_a'] = self.read_no_file(state, mos, no_file, spin=1)
                try:
                    state['sden_b'] = self.read_no_file(state, mos, no_file, spin=-1)
                except IndexError:
                    print("  WARNING: Could not find beta orbitals.")
                    print("  Setting beta=alpha")
                    state['sden_b'] = state['sden_a']
                state['sden']   = state['sden_a'] + state['sden_b']

        for istate, state in enumerate(state_list):
            state['exc_en'] = float(istate + 1) # set fake excitation energy
            state['state_num'] = istate + 1

            state['fname'] = self.ioptions['ana_files'][istate]

            # perform some simplifications to the filename
            # for Columbus:
            tmp_name = state['fname'].replace('MOLDEN/','').replace('./','').replace('molden_no_','').rstrip('.sp')
            tmp_name = tmp_name.replace('state', 'S').replace('drt', 'D')
            # for Q-Chem
            tmp_name = tmp_name.replace('NOs/','').replace('.mo', '')

            tmp_name = os.path.basename(tmp_name)

            state['name'] = tmp_name

        return state_list

    def read_no_file(self, state, ref_mos, no_file, spin=0):
        """
        Read information from a secondary NO file.
        """
        nos = lib_mo.MO_set_molden(file=no_file)
        nos.read(spin=spin)
        #nos.compute_inverse()
        if self.ioptions['rd_ene']:
            nos.set_ens_occs()
        if not self.ioptions['occ_fac'] == 1:
            nos.occs = [occ * self.ioptions['occ_fac'] for occ in nos.occs]
        if self.ioptions['lvprt'] >= 2:
            print("Number of electrons:", sum(nos.occs))

        if not self.ioptions['ignore_irreps'] == []:
            print("  Ignoring irreps: ", self.ioptions['ignore_irreps'])
            for imo, sym in enumerate(nos.syms):
                for iirr in self.ioptions['ignore_irreps']:
                    if iirr in sym:
                        #print "%s -> 0"%sym
                        nos.occs[imo] = 0.

        if not(self.ioptions['min_bf'] == ()):
            ntake = 0
            ibf, vbf = self.ioptions['min_bf']
            print("Selecting only NOs with %s > %f ..."%(nos.bf_labels[ibf], vbf))
            for imo in range(nos.ret_num_mo()):
                pop = nos.ret_mo_pop(imo, dosum=2)
                if pop[ibf] > vbf:
                    if self.ioptions['lvprt'] >= 2: print(("Selecting MO %i, occ = %.3f"%(imo+1, nos.occs[imo])))
                    ntake += 1
                else:
                    nos.occs[imo] = 0.
                    if self.ioptions['lvprt'] >= 2: print(("Discarding MO %i"%(imo+1)))
            print("%i NOs selected."%ntake)

        T = numpy.dot(ref_mos.ret_mo_mat(trnsp=False, inv=True), nos.mo_mat)

        if not self.ioptions['unrestricted']:
            nu_list = [min(occ, 2.-occ) for occ in nos.occs]
            state['nu'] = sum(nu_list)
            state['nu_den'] = numpy.dot(T,
                                  numpy.dot(numpy.diag(nu_list), T.transpose()))

            nunl_list = [occ*occ*(2-occ)*(2-occ) for occ in nos.occs]
            state['nunl'] = sum(nunl_list)
            state['nunl_den'] = numpy.dot(T,
                                  numpy.dot(numpy.diag(nunl_list), T.transpose()))

            nel = sum(nos.occs)
            iy0 = int(nel/2 + 0.5)
            iy1 = int(nel/2 + 0.5) + 1
            try:
                state['y0'] = nos.occs[iy0]
            except IndexError:
                state['y0'] = 0.
            try:
                state['y1'] = nos.occs[iy1]
            except IndexError:
                state['y1'] = 0.

        return numpy.dot(T, numpy.dot(numpy.diag(nos.occs), T.transpose()))

class file_parser_rassi(file_parser_libwfa):
    def read(self, mos):

        if not (self.ioptions['read_libwfa']):
            raise error_handler.MsgError("Please run &WFA within OpenMolcas")

        (energies, oscs, state_list) = self.read_rassi_output(self.ioptions['rfile'])

        return state_list

    def read_rassi_output(self, filen):
        """
        Read the standard output of RASSI.
        """
        state_list = []
        energies = []
        oscs = {}
        libwfa = False

        # This is a bit of a hack to deal with different file names
        ctnumtypes = ['atomic', 'lowdin', 'mulliken']

        rfile = open(filen,'r')

        while(True):
            try:
                line = next(rfile)
            except:
                break

            # if 'Total energies (spin-free)' in line:
            #     words = rfile.next().split()
            #     while(len(words) > 0):
            #         energies.append(float(words[-1]))
            #         words = rfile.next().split()

            if 'SPIN-FREE ENERGIES' in line:
                line = next(rfile)
                if 'Shifted by' in line: next(rfile)
                next(rfile)
                next(rfile)
                words = next(rfile).split()
                while(len(words) > 0):
                    energies.append(float(words[1]))
                    words = next(rfile).split()

            elif '++ Dipole transition strengths' in line:
                next(rfile)
                next(rfile)
                next(rfile)
                next(rfile)
                next(rfile)
                words = next(rfile).split()
                while(len(words) > 1):
                    if words[1] == 'Max':
                        break
                    ist = int(words[0])
                    jst = int(words[1])
                    osc = float(words[2])
                    oscs[(ist, jst)] = osc
                    words = next(rfile).split()
                if not self.ioptions['read_libwfa']:
                    break

            elif 'RASSI analysis for reference state' in line or 'RASSI analysis for state' in line:
                state_list.append({})
                state = state_list[-1]
                istate = len(state_list)

                state['name'] = line.split()[-1]

                line = next(rfile)
                line = next(rfile)
                state['exc_en'] = float(line.split()[1])
                try:
                    state['osc_str'] = oscs[(1, istate)]
                except:
                    pass

                libwfa = True

            elif 'RASSI analysis for transiton from state' in line or 'RASSI analysis for transition from state' in line: # typo fixed
                words = line.split()
                ttyp = words[-1][1:-1]

                # This is a bit of a hack to find the correct file name
                while(len(ctnumtypes)>0):
                    (typ, exctmp, osc, num_at, num_at1, om_at) = self.rmatfile("%s_ctnum_%s.om"%(ttyp, ctnumtypes[0]))
                    if not exctmp == None:
                        if abs(exctmp * units.energy['eV'] - state['exc_en']) > 1.e-3:
                            print(" WARNING: inconsistent energies for %s: %.5f/%.5f - skipping."%(typ, exctmp * units.energy['eV'], state['exc_en']))
                        else:
                            state['Om'] = om_at.sum()
                            state['OmAt'] = om_at
                        break
                    ctnumtypes.pop(0)


            elif libwfa:
                self.parse_line(state, line, rfile)

        rfile.close()

        return energies, oscs, state_list

class file_parser_onetep(file_parser_base):
    """
    Read ONETEP job.
    """
    def read(self, mos):
        self.ioptions['jmol_orbitals'] = False
        self.ioptions['molden_orbitals'] = False

        state_list = []

        pre = self.ioptions['rfile']

        rfile = open("%s.onetep"%pre, 'r')
        while True:
            try:
                line = next(rfile)
            except StopIteration:
                print("Finished parsing %s."%rfile.name)
                break

            if 'Energy (in Ha)' in line:
                while True:
                    line = next(rfile)
                    words = line.split()
                    # print(words[0], words[0] == "Writing")
                    if len(words) != 4 or words[0] == "Writing":
                        break

                    state_list.append({})
                    state = state_list[-1]
                    state['state_ind'] = int(words[0])
                    state['name'] = 'A' + words[0]
                    state['exc_en'] = float(words[1]) * units.energy['eV']
                    state['osc_str'] = float(words[2])

        for istate in range(len(state_list)):
            state = state_list[istate]

            Dao = numpy.zeros([mos.ret_num_bas() // 2, mos.ret_num_bas()], float)
            dmfile = "%s_response_denskern_%i.dkn_dens.mat"%(pre, istate+1)
            for i, line in enumerate(open(dmfile, 'r')):
                words = line.split()
                # Read only the real parts and assert that imaginary parts are zero
                for j in range(0, mos.ret_num_bas(), 2):
                    Dao[j // 2, i] = float(words[j])
                    assert(float(words[j+1]) == 0.)

            state['tden'] = Dao

        return state_list

class file_parser_orca(file_parser_base):
    """
    Read Orca job without using cclib.
    """
    def read(self, mos):
        return self.tden_orca(mos, True)

    def tden_orca(self, mos, rect_dens, filen='orca.cis'):
        """
        Read binary CI vector file from ORCA.
        Authors: S. Mai, F. Plasser
        TODO: delete this part
        """
        print("Reading CI vectors from binary ORCA file %s"%filen)

        CCfile=open(filen,'rb')
        # header
        # consists of 9 4-byte integers, the first 5 of which give useful info
        nvec  =struct.unpack('i', CCfile.read(4))[0]
        print("Number of vectors:", nvec)
        header=[ struct.unpack('i', CCfile.read(4))[0] for i in range(8) ]
        #print 'header:', header

        # header array contains:
        # [0] index of first alpha occ,  is equal to number of frozen alphas
        # [1] index of last  alpha occ
        # [2] index of first alpha virt
        # [3] index of last  alpha virt, header[3]+1 is equal to number of bfs
        # [4] index of first beta  occ,  for restricted equal to -1
        # [5] index of last  beta  occ,  for restricted equal to -1
        # [6] index of first beta  virt, for restricted equal to -1
        # [7] index of last  beta  virt, for restricted equal to -1

        if any( [ header[i]!=-1 for i in range(4,8) ] ):
            #raise error_handler.MsgError("No support for unrestricted MOs")
            print("Detected unrestricted calculation!")
            restr=False
        else:
            restr=True

        if restr:
            do_Alpha = True
        else:
            do_Alpha = (self.ioptions['spin'] >= 0)
        print("! Reading %s coefficients ..." % (["BETA","ALPHA"][do_Alpha]))

        NFA = header[0]
        NOA = header[1] - header[0] + 1
        NVA = header[3] - header[2] + 1
        NFB = header[4]
        NOB = header[5] - header[4] + 1
        NVB = header[7] - header[6] + 1

        if do_Alpha:
          nfrzc = NFA
          nocc  = NFA + NOA
          nact  = NOA
          nvir  = NVA
          nmo   = NFA + NOA + NVA
          lenci = NOA * NVA
          if restr:
            skipci = 0
          else:
            skipci = NOB * NVB
        else:
          nfrzc = NFB
          nocc  = NFB + NOB
          nact  = NOB
          nvir  = NVB
          nmo   = NFB + NOB + NVB
          lenci = NOB * NVB
          skipci = NOA * NVA

        print('  nmo: %i , nocc: %i , nact: %i , nvir: %i'%(nmo,nocc,nact,nvir))
        if nmo != mos.ret_num_mo():
            raise error_handler.MsgError("Inconsistent number of MOs")

        # loop over states
        # for non-TDA order is: X+Y of 1, X-Y of 1, X+Y of 2, X-Y of 2, ...
        # triplets come after singlets, observe the multiplicity!
        # for unrestricted: each vector has first alpha part, then beta part
        # so if we only analyze one spin, then we have to skip the other
        rootinfo=[]
        istate = -1
        irootdict = {}
        TDA=True
        energies = [] # In the case of RPA, the energies are written with the first N entries
        state_list = []
        for ivec in range(nvec):
            # header of each vector
            # contains 6 4-byte ints, then 1 8-byte double, then 8 byte unknown
            d0,d1,mult,d2,iroot,d3 = struct.unpack('iiiiii', CCfile.read(24))
            rootinfo.append( (mult,iroot) )
            ene,d3 = struct.unpack('dd', CCfile.read(16))
            print('  mult: %i , iroot: %i'%(mult,iroot))
            #print(d1,d2,d3)
            print('  Ene: %.4f eV'%(ene*units.energy['eV']))
            energies.append(ene)
            # -> energy looks ok but the order is strange for RPA

            # then comes nact * nvirt 8-byte doubles with the coefficients
            if not restr and not do_Alpha:
              CCfile.read(skipci*8)
            coeff = struct.unpack(lenci*'d', CCfile.read(lenci*8))
            if not restr and do_Alpha:
              CCfile.read(skipci*8)

            if ivec==1 and (mult,iroot)==rootinfo[0]:
              TDA=False
              print("Detected a non-TDA calculation!")

            if ivec>=1 and mult!=rootinfo[-2][0]:
              triplets=True
              #istate=-1


            if TDA or ivec%2==0:
                istate+=1
                irootdict[mult] = irootdict.get(mult, -1) + 1
                state_list.append({})
                state = state_list[-1]
                state['state_ind'] = irootdict[mult] + 1
                state['mult'] = mult
                state['irrep'] = 'A'
                state['name'] = '%i(%i)%s'%(state['state_ind'], state['mult'] ,state['irrep'])
                state['exc_en'] = energies.pop(0)*units.energy['eV']
                state['tden'] = self.init_den(mos, rect=rect_dens)

                state['tden'][nfrzc:nocc, nocc:nmo] = numpy.reshape(coeff, [nact,nvir])
            else:
            # in this case, we have a non-TDA state!
            # and we need to compute (prevvector+currentvector)/2 = X vector
            # this is added to the vector from the previous "if"
                print('Constructing X-vector of RPA state')
                state['tden'][nfrzc:nocc, nocc:nmo] += numpy.reshape(coeff, [nact,nvir])
                state['tden'][nfrzc:nocc, nocc:nmo] *= .5

        return state_list

class file_parser_orca_CAS(file_parser_base):

    def read(self, mos):
        return self.tden_orca(mos, False)

    def tden_orca(self, mos, rect_dens, filen1='orca.densities', filen2='orca.densitiesinfo'):
        '''
        as far as i found out by looking into multiple files produced by orca6.1 each section has length of 661 bytes starting
        with a nul terminated string representing that densities name. The header seed to variied in length. It might be also due 
        the presents of a nul terminated string starting at byte 4 in the data.
        The first 4 bytes seem to be a uint32/int32 stating the number of available densities.

        '''
        state_list = []

        outFile = open(self.ioptions['rfile'], 'r')
        startRead = False
        for l in outFile:
            if startRead:
                if l =="\n":
                    startRead = False
                    continue

                words = l.split()
                state = {'mult': int(words[2].split("-")[1][0]), 'irrep': words[2].split("-")[1][1:], 'state_ind': int(words[2].split("-")[0]), 'exc_en': float(words[3])}
                state['name'] = '%i(%i)%s'%(state['state_ind'], state['mult'] ,state['irrep'])
                state_list.append(state)
                
            elif "                     ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS    \n" == l:
                [next(outFile) for i in range(4)]
                state_list = []
                startRead=True
            

        dens_info = open(filen2, 'rb')
        data = dens_info.read()
        dens_info.close()
        file_size = len(data)
        n_sections = struct.unpack('<i', data[:4])[0]
        header_len = file_size - 661 * n_sections
        header = data[:header_len]

        parse_instructions = []
        for i, d in enumerate([data[i : i + 661] for i in range(header_len, len(data), 661)]):
            name_end = d.find(b'\x00')
            
            if name_end != -1:
                section_name = d[0:name_end].decode('utf-8', errors='ignore')
            else:
                section_name = "UNKNOWN_SECTION"

            descript_pos = d.find(b'\xff\xff\xff\xff\xff\xff\xff\xff') + 8

            ambiguous_bytes = d[descript_pos + 8 : descript_pos + 13]     
            if ambiguous_bytes != b'\x08\x00\x00\x00\x00':
                print(f'Notice: Ambiguous bytes changed at section {i}: {ambiguous_bytes.hex()}')
            
            # in all files i looked the 661 byte long regions starts a 8 byte long FF terminater starting at byte 520, than based on the nuls in the hex
            # there seed to be multiple int32s and int64s in that region afterwords. The only important here are the first two int32 after the terminator
            # as these provide the size of the corresponding matrix in .densities. i am not sure which is row and which is call as in all my files they
            # where square matrixes, i still parse out all the numbers i think i found but i will only use the two mentioned
            int32_block = d[descript_pos : descript_pos + 8] + d[descript_pos + 13 : descript_pos + 37]
            int64_block = d[descript_pos + 37:]
            assert len(int32_block) == 32, len(int32_block)
            assert len(int64_block) == 96, len(int64_block)
            
            unpacked_uint32s = struct.unpack('<IIIIIIII', int32_block)
            unpacked_uint64s = struct.unpack('<QQQQQQQQQQQQ', int64_block)
            
            # selecting only the transition densities starting from the ground root 
            if ("Tdens-CAS-" in section_name) and (section_name.split('-')[4]=="0") and (section_name.split('-')[5]!="0"):
                parse_instructions.append((unpacked_uint32s[0], unpacked_uint32s[1]))
                print(section_name)
            else:
                parse_instructions.append(unpacked_uint32s[0] * unpacked_uint32s[1])


        # starting here i will actually parse the .densities file. When using a double for each value expected by the matrix sizes gained from the the
        # densitiesinfo, the length of the densities file fits exactly.
        bytes_per_value = 8
        densfile = open(filen1, 'rb')

        print(f"found {len(state_list)} states")
        i = 0
        for item in parse_instructions:
            if isinstance(item, int):
                densfile.seek(bytes_per_value * item, 1)
            elif isinstance(item, tuple):
                rows, cols = item # might be flipped. In my data only square matrixes where present."
                bytes_to_read = bytes_per_value * rows * cols 
                raw = densfile.read(bytes_to_read)

                if len(raw) != bytes_to_read:
                    raise EOFError(f"Unexpected end of file while reading matrix of size {rows}x{cols}. Expected {bytes_to_read} bytes, but got {len(raw)} bytes.")

                D_AO = numpy.frombuffer(raw, dtype=numpy.float64).reshape((rows,cols))
                temp   = mos.CdotD(D_AO.T, trnsp=False, inv=True)   # C⁻¹ @ D_AOᵀ      → (nmo × nbas)
                D_MO   = mos.MdotC(temp, trnsp=True, inv=True)      # C⁻¹ @ D_AOᵀ @ C⁻ᵀ  → (nmo × nmo)
                state_list[i]['tden'] = D_MO
                i += 1

        # just doing some sanity check that everything was as expected and we now reached the end of the densities file.
        extra_byte = densfile.read(1)
        if extra_byte != b'':
            current_pos = densfile.tell() - 1  
            densfile.seek(0, 2)               
            total_size = densfile.tell()
            remaining_bytes = total_size - current_pos
            
            raise ValueError(f"Layout fully processed, but the end of the file was not reached! There are still {remaining_bytes} unparsed bytes left in '{filen1}'.")
            
        return state_list
