"""
This is an interface to the external cclib Library.
http://cclib.github.io
Download and install cclib if you want to use the functions.
"""

from __future__ import print_function, division
import struct
import numpy
from . import file_parser, lib_mo, error_handler, units, lib_struc
try:
    from openbabel import openbabel
except ImportError:
    print(" *** Warning: python-openbabel not found! ***")
    print(" Using emulation program with limited capabilities ...")
    from . import OB_repl as openbabel

class file_parser_cclib(file_parser.file_parser_base):
    def __init__(self, *args, **kwargs):
        file_parser.file_parser_base.__init__(self, *args, **kwargs)

        self.data, self.prog = self.get_data(self.ioptions['rtype'], self.ioptions['rfile'])

    def get_data(self, rtype, rfile):
        """
        Call cclib to read in the data.
        """
        try:
            import cclib.parser
            print(("Using cclib version %s"%cclib.__version__))
        except ImportError:
            print("\n ERROR: Did not find the external cclib package!")
            print("Please, download and install cclib to use this functionality.\n")
            raise

        cfile = cclib.parser.ccopen(rfile)
        if cfile == None:
            raise error_handler.MsgError('File %s cannot be parsed by cclib!'%rfile)
        else:
            print(("Parsing %s ..."%str(cfile)))

        prog = str(cfile).split()[0]

        return cfile.parse(), prog

    def read(self, mos, rect_dens=True):
        """
        Convert the information to data analyzed by TheoDORE.
        """
        state_list = []

        try:
            core_orbs = sum(self.data.coreelectrons) / 2
            if core_orbs > 0:
                print(("\n WARNING: %i core orbitals detected"%core_orbs))
                print(" Please check the results carefully!\n")
                rect_dens = False
        except AttributeError:
            pass

        for ist, exc_en in enumerate(self.data.etenergies):
            state_list.append({})
            state = state_list[-1]

            state['exc_en'] = exc_en / units.energy['rcm'] * units.energy['eV']
            try:
                state['osc_str'] = self.data.etoscs[ist]
            except AttributeError:
                state['osc_str'] = -1.
            try:
                state['irrep'] = self.data.etsyms[ist].replace('let-', '').replace('Not specified', 'X')
            except AttributeError:
                state['irrep'] = 'X'

        if self.prog == 'ADF':
            self.tden_adf(state_list, mos, rect_dens)
        elif self.prog == 'ORCA':
            if self.ioptions['read_binary']:
                self.tden_orca(state_list, mos, rect_dens)
            else:
                print('Reading ORCA input from stdout (read_binary=False)')
                self.tden_cclib(state_list, mos, rect_dens)
        else:
            self.tden_cclib(state_list, mos, rect_dens)

        self.check_RKS()

        return state_list

    def tden_cclib(self, state_list, mos, rect_dens):
        for ist, state in enumerate(state_list):
            state['state_ind'] = ist + 1
            state['name'] = '%i%s'%(state['state_ind'],state['irrep'])
            print(("\n" + state['name']))

            state['tden'] = self.init_den(mos, rect=rect_dens)
            for conf in self.data.etsecs[ist]:
                [(iocc, spocc), (ivirt, spvirt), coeff] = conf
                if self.ioptions['spin'] >= 0:
                    if spocc == 1:
                        continue
                if self.ioptions['spin'] == -1:
                    if spocc == 0:
                        continue
                assert(iocc < ivirt) # de-excitations cannot be handled at this point
                state['tden'][iocc, ivirt] = coeff

                if coeff*coeff > 0.05:
                    print((" (%i->%i), coeff=% .4f"%(iocc+1, ivirt+1, coeff)))
                elif numpy.isnan(coeff):
                    print((" ERROR: (%i->%i), coeff=% .4f"%(iocc+1, ivirt+1, coeff)))
                    if self.prog == 'ORCA':
                        print("For ORCA/RPA jobs, please set read_binary=True\n")
                    raise error_handler.MsgError("Not-a-number encountered")

    def tden_orca(self, state_list, mos, rect_dens, filen='orca.cis'):
        """
        Read binary CI vector file from ORCA.
        Authors: S. Mai, F. Plasser
        TODO: delte this subroutine (it is moved to file_parser.py)
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
        iroot2 = -1
        TDA=True
        for ivec in range(nvec):
            # header of each vector
            # contains 6 4-byte ints, then 1 8-byte double, then 8 byte unknown
            d0,d1,mult,d2,iroot,d3 = struct.unpack('iiiiii', CCfile.read(24))
            rootinfo.append( (mult,iroot) )
            ene,d3 = struct.unpack('dd', CCfile.read(16))
            print('  mult: %i , iroot: %i'%(mult,iroot))
            #print d1,d2,d3
            #print '  Ene: %.4f, Osc: %.4f'%(ene*units.energy['eV'], osc)
            # -> energy does not look consistent

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
              iroot2=-1


            if TDA or ivec%2==0:
                istate+=1
                iroot2+=1
            #if prevroot!=iroot:
                state = state_list[istate]
                state['state_ind'] = iroot2 + 1
                state['mult'] = mult
                state['irrep'] = 'A'
                state['name'] = '%i(%i)%s'%(state['state_ind'], state['mult'] ,state['irrep'])
                state['tden'] = self.init_den(mos, rect=rect_dens)
                #assert abs(state['exc_en']-ene*units.energy['eV']) < 1.e-4, 'Incorrect energies'
                # -> does not work for all states

                state['tden'][nfrzc:nocc, nocc:nmo] = numpy.reshape(coeff, [nact,nvir])
            else:
            # in this case, we have a non-TDA state!
            # and we need to compute (prevvector+currentvector)/2 = X vector
                print('Constructing X-vector of RPA state')
                state['tden'][nfrzc:nocc, nocc:nmo] += numpy.reshape(coeff, [nact,nvir])
                state['tden'][nfrzc:nocc, nocc:nmo] *= .5

    def check(self, lvprt=1, maxerr=50):
        """
        Check if the input file can be used by cclib.
        lvprt  ... print level
        maxerr ... maximum allowed error code
        """
        errcode = 0

        if lvprt >= 1:
            print(("\nChecking if the logfile %s can be parsed by cclib ..."%self.ioptions['rfile']))
            print("\nEssential attributes:")

        for attr in ['mocoeffs', 'atombasis', 'natom', 'homos', 'moenergies', 'etenergies', 'etsyms', 'etsecs']:
            chk = hasattr(self.data, attr)
            if not chk: errcode = 3

            if lvprt >= 1:
                print(('%15s ... %s'%(attr, chk)))
            if chk and lvprt >= 2:
                print(getattr(self.data, attr))

        if lvprt >= 1:
            print("\nOptional attributes:")

        for attr in ['etoscs', 'aooverlaps', 'mosyms']:
            chk = hasattr(self.data, attr)

            if lvprt >= 1:
                print(('%15s ... %s'%(attr, chk)))
            if chk and lvprt >= 2:
                print(getattr(self.data, attr))

        if lvprt >= 1:
            print("\nAttributes for structure parsing and creation of Molden file:")

        if self.prog.lower() in ['orca']:
            print(" Conversion to Molden format not supported for %s!"%self.prog)
            errcode = max(1, errcode)

        for attr in ['gbasis', 'natom', 'atomcoords', 'atomnos']:
            chk = hasattr(self.data, attr)
            if not chk: errcode = max(1, errcode)

            if lvprt >= 1:
                print(('%15s ... %s'%(attr, chk)))
            if chk and lvprt >= 2:
                print(getattr(self.data, attr))

        if errcode > maxerr:
            raise error_handler.MsgError("The file cannot be parsed by cclib")

        return errcode

    def check_RKS(self, lvprt=1):
        """
        Check if this is an RKS job.
        """
        chk = (len(self.data.homos) == 1)

        if chk:
            if lvprt >= 1:
                print("\n Detected RHF/RKS calculation.")
            if self.ioptions['spin'] != 0:
                raise error_handler.MsgError("RHF/RKS calculation")
        else:
            if lvprt >= 1:
                print("\n Using UHF/UKS mode.")
            if self.ioptions['spin'] == 0:
                raise error_handler.MsgError("Use analyze_tden_unr.py for unrestricted calculations.")

        return chk

    def read_mos(self):
        """
        Return an MO_set object in TheoDORE format.
        """
        if self.prog == 'ADF':
            #mos = MO_set_adf(file = None)
            mos = MO_set_cclib(file = None)
        else:
            mos = MO_set_cclib(file = None)

        do_Alpha = (self.ioptions['spin'] >= 0)
        mos.read(self.data, lvprt=2, do_Alpha=do_Alpha)

        return mos

    def read_mos_from(self, mo_rfile, lvprt=1, en_thresh=1.e-4, coor_thresh=1.e-4):
        """
        Read the MOs from a separate log file (parsed with cclib), e.g. for
        Gaussian TDDFT jobs without pop=full where the orbitals were read
        (guess=read) from a single point with pop=full.
        The two files are checked for consistency: same atoms and geometry,
        same SCF energy and number of basis functions, and the same occupied
        orbital energies if both files contain them.
        """
        print("\nReading MOs from separate file %s"%mo_rfile)
        mo_data, mo_prog = self.get_data(self.ioptions['rtype'], mo_rfile)

        if not hasattr(mo_data, 'mocoeffs'):
            raise error_handler.MsgError("No MO coefficients in %s (use pop=full for Gaussian)"%mo_rfile)

        errors = []
        if hasattr(self.data, 'nbasis') and self.data.nbasis != mo_data.nbasis:
            errors.append("number of basis functions: %i / %i"%(self.data.nbasis, mo_data.nbasis))
        if not numpy.array_equal(self.data.atomnos, mo_data.atomnos):
            errors.append("atoms differ")
        else:
            dcoor = numpy.max(abs(numpy.array(self.data.atomcoords[-1]) - numpy.array(mo_data.atomcoords[-1])))
            if dcoor > coor_thresh:
                errors.append("geometries differ by up to %.2e A (same orientation needed)"%dcoor)
        if hasattr(self.data, 'scfenergies') and hasattr(mo_data, 'scfenergies'):
            # scfenergies are in eV in cclib
            den = abs(self.data.scfenergies[-1] - mo_data.scfenergies[-1]) / units.energy['eV']
            if den > en_thresh:
                errors.append("SCF energies differ by %.2e Eh"%den)
        if hasattr(self.data, 'moenergies') and hasattr(mo_data, 'moenergies'):
            nocc = min(self.data.homos[0], mo_data.homos[0]) + 1
            if self.data.homos[0] != mo_data.homos[0]:
                errors.append("HOMO index: %i / %i"%(self.data.homos[0], mo_data.homos[0]))
            elif nocc > 0:
                den = numpy.max(abs(numpy.array(self.data.moenergies[0][:nocc]) - numpy.array(mo_data.moenergies[0][:nocc]))) / units.energy['eV']
                if den > en_thresh:
                    errors.append("occupied orbital energies differ by up to %.2e Eh"%den)

        if len(errors) > 0:
            raise error_handler.MsgError("%s and %s are not consistent:\n  "%(self.ioptions['rfile'], mo_rfile) + "\n  ".join(errors))
        if lvprt >= 1:
            print(" Consistency of %s and %s checked (geometry, SCF energy, orbital energies)"%(self.ioptions['rfile'], mo_rfile))

        mos = MO_set_cclib(file = None)
        do_Alpha = (self.ioptions['spin'] >= 0)
        mos.read(mo_data, lvprt=2, do_Alpha=do_Alpha)

        return mos

    def ret_struc(self,lvprt=1):
        struc = structure_cclib()
        struc.read_cclib(self.data)

        return struc

class MO_set_cclib(lib_mo.MO_set_molden):
    def read(self, data, lvprt=1, do_Alpha=True):
        """
        Read cclib data and convert to TheoDORE format.
        """

        if do_Alpha:
          self.mo_mat = data.mocoeffs[0].transpose()
          print("! Reading %s orbitals ..." % (["BETA","ALPHA"][do_Alpha]))
        else:
          if len(data.mocoeffs)>1:
            self.mo_mat = data.mocoeffs[1].transpose()
            print("! Reading %s orbitals ..." % (["BETA","ALPHA"][do_Alpha]))
          else:
            raise error_handler.MsgError("No beta MO-coefficients available")
        if lvprt >= 1:
            print(("MO-matrix read from cclib, dimension: %i x %i"%(len(self.mo_mat), len(self.mo_mat[0]))))

        self.num_at = data.natom

        self.ens = data.moenergies[0]
        self.ihomo = data.homos[0]
        self.occs = (self.ihomo + 1) * [1.] + (len(self.ens) - self.ihomo - 1) * [0.]

        self.basis_fcts = [lib_mo.basis_fct(-1) for ibas in range(self.ret_num_bas())]
        maxbas = -1
        for iat, baslist in enumerate(data.atombasis):
            maxbas = max([maxbas] + baslist)
            for ibas in baslist:
                self.basis_fcts[ibas].set(iat + 1)

        assert(self.ret_num_bas() == maxbas+1)

        try:
            self.syms = data.mosyms
        except AttributeError:
            self.syms = ['X'] * len(self.ens)

        try: self.S = data.aooverlaps
        except AttributeError: pass

        try: self.atomnos = data.atomnos
        except AttributeError: pass

        try: self.coords = data.atomcoords[-1]
        except AttributeError: pass

        try: self.gbasis = data.gbasis
        except AttributeError: pass

    def ret_ihomo(self):
        return self.ihomo

    def write_molden_file(self, fname='out.mld', cfmt='% 10E', occmin=-1):
        """
        Write a file in Molden format.
        """
        self.header = "[Molden Format]\n"

        self.header+= "[Atoms] Angs\n"
        for iat in range(self.num_at):
            atno = self.atomnos[iat]
            (x, y, z) = self.coords[iat]
            self.header+= "%3s %5i %3i %12.6f %12.6f %12.6f\n"%\
                (lib_struc.Z_symbol_dict[atno], iat+1, atno, x, y, z)

        # Count the basis functions to see if it spherical or Cartesian
        numl_s={'s':1,'p':3,'sp':4,'d':5,'f': 7,'g': 9}
        numl_c={'s':1,'p':3,'sp':4,'d':6,'f':10,'g':15}
        nbas = [0,0]
        self.header+= "[GTO]\n"
        for iat, at in enumerate(self.gbasis):
            self.header+= "%3i 0\n"%(iat+1)
            for bas in at:
                typ = bas[0]
                exps = bas[1]
                nbas[0] += numl_s[typ.lower()]
                nbas[1] += numl_c[typ.lower()]

                self.header+= "%3s %5i 1.00\n"%(typ, len(exps))

                for exp in exps:
                    self.header+= "%18.10E %18.10E\n"%exp
            self.header+= "\n"

        if nbas[1] == self.ret_num_bas():
            print("Assuming Cartesian basis functions")
        elif nbas[0] == self.ret_num_bas():
            print("Assuming spherical basis functions")
            self.header += '\n[5D]\n[7F]\n[9G]\n'
        else:
            print("\n WARNING: Inconsistent number of basis functions!")
            print("Spherical: %i, Cartesian: %i, actual: %i"%(nbas[0], nbas[1], self.ret_num_bas()))

        self.export_AO(self.ens, self.occs, self.mo_mat.transpose(), fname, cfmt, occmin)

class structure_cclib(lib_struc.structure):
    def read_cclib(self, data, ind=-1, lvprt=1):
        if lvprt>=1:
            print(("Reading cclib structure with %i atoms."%data.natom))

        self.mol = openbabel.OBMol()

        for iat in range(data.natom):
            obatom = openbabel.OBAtom()
            obatom.SetAtomicNum(int(data.atomnos[iat]))
            coords = data.atomcoords[ind][iat]
            obatom.SetVector(*coords)

            self.mol.AddAtom(obatom)
