"""
Interface to the JSON files written by orca_2json (ORCA 6).

Run the excited-state job and convert the gbw file:
    orca_2json job
which writes job.json with the MO coefficients, the basis function labels and
the CIS/TDA/TDDFT amplitudes ("Xia"). If job.property.json exists next to
job.json, oscillator strengths are read from it as well.

No external libraries are needed, and the binary .cis file is not used, so the
interface does not depend on the binary file layout.
"""

from __future__ import print_function, division
import json
import os
import numpy
from . import file_parser, lib_mo, error_handler, units


def load_orca_json(fname):
    """
    Read a JSON file written by orca_2json and check that it has a molecule.
    """
    with open(fname, 'r') as f:
        data = json.load(f)
    if not 'Molecule' in data:
        raise error_handler.MsgError("%s does not contain a 'Molecule' section (written by orca_2json?)"%fname)
    return data

def orca_version(data):
    try:
        return data['ORCA Header']['Version']
    except (KeyError, TypeError):
        return 'unknown'

def is_unrestricted(mol):
    return mol.get('HFTyp', 'RHF').upper() in ['UHF', 'UKS']

class MO_set_orca_json(lib_mo.MO_set):
    """
    MO coefficients, energies and basis function to atom map from orca_2json.
    """
    def read(self, data=None, lvprt=1, spin=0):
        """
        data - dictionary from load_orca_json (read from self.file if None)
        spin - 0 restricted, 1 alpha, -1 beta (unrestricted only)
        """
        if data is None:
            data = load_orca_json(self.file)
        mol = data['Molecule']

        if lvprt >= 1:
            print("Reading MOs from orca_2json file (ORCA version %s)"%orca_version(data))

        morb = mol['MolecularOrbitals']
        mos = morb['MOs']
        labels = morb['OrbitalLabels']
        num_bas = len(labels)

        if is_unrestricted(mol):
            if len(mos) % 2 != 0:
                raise error_handler.MsgError("Odd number of MOs for an unrestricted calculation")
            nmo = len(mos) // 2
            if spin == -1:
                mos = mos[nmo:]
                if lvprt >= 1: print("! Reading BETA orbitals ...")
            else:
                mos = mos[:nmo]
                if lvprt >= 1: print("! Reading ALPHA orbitals ...")
        elif spin == -1:
            raise error_handler.MsgError("No beta MO-coefficients available (restricted calculation)")

        self.mo_mat = numpy.array([mo['MOCoefficients'] for mo in mos], float).T
        if self.mo_mat.shape[0] != num_bas:
            raise error_handler.MsgError("Inconsistent number of basis functions: %i labels, %i coefficients"%(num_bas, self.mo_mat.shape[0]))
        self.ens  = [mo['OrbitalEnergy'] for mo in mos]
        self.occs = [mo['Occupancy'] for mo in mos]
        self.syms = [mo.get('OrbitalSymLabel', 'X') for mo in mos]

        # basis function labels, e.g. '0O   1s', '1H   2pz' -> atom index, l, ml
        self.basis_fcts = []
        self.bf_labels = []
        for lab in labels:
            words = lab.split()
            iat = int(''.join(c for c in words[0] if c.isdigit()))
            shell = words[1].lstrip('0123456789')
            self.basis_fcts.append(lib_mo.basis_fct(iat + 1, shell[:1], shell[1:]))
            if not self.basis_fcts[-1].label() in self.bf_labels:
                self.bf_labels.append(self.basis_fcts[-1].label())

        # atoms (Angstrom)
        unit = mol.get('CoordinateUnits', 'Angs')
        fac = units.length['A'] if unit.lower().startswith(('b', 'au')) else 1.
        self.at_dicts = []
        for atom in mol['Atoms']:
            (x, y, z) = atom['Coords']
            self.at_dicts.append({'Z': int(atom['ElementNumber']), 'x': x*fac, 'y': y*fac, 'z': z*fac})
        self.num_at = len(self.at_dicts)

        # Overlap matrix, only present if requested in the orca_2json configuration
        if 'S-Matrix' in mol:
            S = numpy.array(mol['S-Matrix'], float)
            if S.shape == (num_bas, num_bas):
                self.S = S

        if lvprt >= 1:
            print("MO-matrix read from JSON, dimension: %i x %i"%self.mo_mat.shape)
            print("Number of atoms: %i"%self.num_at)

class file_parser_orca_json(file_parser.file_parser_base):
    """
    Read the CIS/TDA/TDDFT amplitudes from an orca_2json file.
    The excitation vectors are read from the blocks that contain a list of
    roots ("SingletRoots", "TripletRoots", ...) together with the orbital
    window "OrbWin". For RPA the X amplitudes are used (as for orca.cis).
    """
    mult_dict = {'Singlet': 1, 'Triplet': 3}

    def read(self, mos, rect_dens=True, data=None):
        rfile = self.ioptions['rfile']
        if data is None:
            data = load_orca_json(rfile)
        mol = data['Molecule']

        spin = self.ioptions.get('spin', strict=False)
        if spin is None: spin = 0
        ispin = 1 if spin == -1 else 0

        roots = self.collect_roots(mol)
        if len(roots) == 0:
            raise error_handler.MsgError("No excited-state amplitudes found in %s"%rfile)

        nmo = mos.ret_num_mo()
        state_list = []
        irootdict = {}
        any_Y = False
        for (mult, orbwin, root) in roots:
            win = self.orbital_window(orbwin, ispin)
            (ifrz, ilocc, ifvirt, ilvirt) = win
            nact  = ilocc - ifrz + 1
            nvirt = ilvirt - ifvirt + 1
            if ilvirt >= nmo:
                raise error_handler.MsgError("Orbital window %s does not fit %i MOs"%(str(win), nmo))

            X = self.spin_block(root['Xia'], ispin)
            if X.shape != (nact, nvirt):
                raise error_handler.MsgError("Amplitudes with shape %s do not match the orbital window %s"%(str(X.shape), str(win)))
            if 'Yia' in root:
                any_Y = True

            irootdict[mult] = irootdict.get(mult, 0) + 1
            state = {}
            state['state_ind'] = irootdict[mult]
            state['mult'] = mult
            state['irrep'] = self.irrep_label(mol, root)
            state['name'] = '%i(%i)%s'%(state['state_ind'], state['mult'], state['irrep'])
            state['exc_en'] = root['Energy'] * units.energy['eV']
            state['tden'] = self.init_den(mos, rect=rect_dens)
            if rect_dens and ilocc >= state['tden'].shape[0]:
                raise error_handler.MsgError("Orbital window %s inconsistent with the MO occupations"%str(win))
            state['tden'][ifrz:ilocc+1, ifvirt:ilvirt+1] = X
            state_list.append(state)

        if any_Y:
            print(" Note: RPA calculation, using the X amplitudes (as for orca.cis)")

        self.read_osc(state_list, rfile)

        return state_list

    def collect_roots(self, mol):
        """
        Return a list of (multiplicity, OrbWin, root) for all roots with
        amplitudes, in the order of the file.
        """
        roots = []
        for key, block in mol.items():
            if not isinstance(block, dict) or not 'OrbWin' in block:
                continue
            for rkey, rlist in block.items():
                if not rkey.endswith('Roots') or not isinstance(rlist, list):
                    continue
                mult = self.mult_dict.get(rkey[:-5], mol.get('Multiplicity', 1))
                for root in rlist:
                    if isinstance(root, dict) and 'Xia' in root:
                        roots.append((mult, block['OrbWin'], root))
        return roots

    def orbital_window(self, orbwin, ispin):
        """
        OrbWin is [first occ, last occ, first virt, last virt], given either
        once or once per spin. Indices start at 0.
        """
        if len(orbwin) > 0 and isinstance(orbwin[0], list):
            win = orbwin[ispin]
        else:
            win = orbwin
        win = [int(i) for i in win]
        if len(win) != 4 or min(win) < 0:
            raise error_handler.MsgError("No orbital window for this spin: %s"%str(orbwin))
        return win

    def spin_block(self, amp, ispin):
        """
        The amplitudes are stored as [spin][occ][virt].
        """
        if len(amp) > 0 and isinstance(amp[0], list) and len(amp[0]) > 0 and isinstance(amp[0][0], list):
            if ispin >= len(amp):
                raise error_handler.MsgError("No amplitudes for spin index %i"%ispin)
            return numpy.array(amp[ispin], float)
        return numpy.array(amp, float)

    def irrep_label(self, mol, root):
        irrep = root.get('Irrep', 0)
        if isinstance(irrep, str):
            return irrep
        if mol.get('PointGroup', 'C1').upper() == 'C1':
            return 'A'
        return 'irrep%i'%irrep

    def read_osc(self, state_list, rfile, thresh=1.e-4):
        """
        Read oscillator strengths from the property file (job.property.json),
        matching the states by their excitation energy.
        """
        if rfile.endswith('.json'):
            pfile = rfile[:-5] + '.property.json'
        else:
            return
        if not os.path.exists(pfile):
            return

        with open(pfile, 'r') as f:
            prop = json.load(f)
        spectra = []
        for geom in prop.get('Geometries', []):
            for spec in geom.get('Absorption_Spectrum', []):
                if spec.get('Representation', 'Length') == 'Length':
                    spectra.append(spec)
        if len(spectra) == 0:
            return

        rows = spectra[-1].get('ExcitationEnergies', [])
        mults = spectra[-1].get('Multiplicities', [])
        print("Reading oscillator strengths from %s"%pfile)
        used = set()
        for state in state_list:
            for irow, row in enumerate(rows):
                if irow in used or abs(row[0] - state['exc_en']) >= thresh:
                    continue
                if irow < len(mults) and int(round(mults[irow][-1])) != state['mult']:
                    continue
                state['osc_str'] = row[3]
                used.add(irow)
                break
