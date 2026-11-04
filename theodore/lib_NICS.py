"""
Analysis routines in the context of NICS (nucleus independent chemical shift)
calculations.

Author: Felix Plasser
Date: Nov. 2020
"""

import re
import numpy
from theodore import units

# numpy print settings used for the diagnostic output of NICS_point.diag().
# They are set explicitly because other TheoDORE modules change the global
# numpy print options on import, which made the output order dependent.
NICS_PRINTOPTIONS = {'precision': 8, 'suppress': False, 'linewidth': 75}

# Ring geometry and ghost atom positions
#   These helpers are general and do not depend on the quantum chemistry code.

def ring_frame(ring_coords):
    """
    Return the centroid and the unit normal of a ring.
    ring_coords - (n x 3) array of the ring atom coordinates (n >= 3)
    The normal is the direction of smallest extent of the atoms (least-squares
    plane). Its sign follows the right-hand rule with respect to the order in
    which the ring atoms are given, so the same atom order always gives the
    same side of the ring.
    """
    coords = numpy.array(ring_coords, float)
    if coords.ndim != 2 or coords.shape[1] != 3 or coords.shape[0] < 3:
        raise ValueError("ring_frame needs at least three 3D coordinates")
    centroid = numpy.mean(coords, axis=0)
    rel = coords - centroid
    normal = numpy.linalg.svd(rel)[2][2]

    # Orientation: sum of the cross products of consecutive atoms
    orient = numpy.zeros(3)
    for i in range(len(rel)):
        orient += numpy.cross(rel[i], rel[(i+1) % len(rel)])
    if numpy.dot(orient, normal) < 0.:
        normal = -normal

    return centroid, normal

def ring_probe_points(ring_coords, heights=(0., 1.), both_sides=False):
    """
    Return the positions of ghost atoms on the ring axis.
    heights    - distances (in the units of ring_coords, usually Angstrom)
                 from the ring centroid along the normal, e.g. (0., 1.) for
                 NICS(0) and NICS(1)
    both_sides - also place the points at the negative heights
    """
    centroid, normal = ring_frame(ring_coords)
    hlist = []
    for h in heights:
        hlist.append(float(h))
        if both_sides and h != 0.:
            hlist.append(-float(h))
    return [centroid + h * normal for h in hlist], hlist

def ring_scan_points(ring_coords, hmin=0., hmax=3., step=0.5):
    """
    Return equally spaced ghost atom positions on the ring axis between
    hmin and hmax (inclusive), e.g. for a NICS scan.
    """
    nstep = int(round((hmax - hmin) / step))
    heights = [hmin + i * step for i in range(nstep + 1)]
    return ring_probe_points(ring_coords, heights)

# Ghost atom input (specific to the quantum chemistry code)
#   Only Gaussian is implemented so far.

def gaussian_bq_lines(points, fmt='%12.6f'):
    """
    Return the lines of a Gaussian molecule specification for ghost atoms (Bq).
    """
    # avoid printing -0.000000
    return ['Bq ' + ' '.join(fmt%(x if abs(x) >= 5.e-7 else 0.) for x in point) for point in points]

# Output

class NICS_parser:
    """
    Parse and visualise NICS calculations.
    """
    def __init__(self, logfile=None):
        self.NICS_data = []
        self.atoms = [] # real atoms: [{'Z': , 'coor': numpy.array}], if parsed
        self.mol_ind = 0
        if not logfile is None:
            self.read(logfile)

    def read(self, logfile):
        """
        Read data from quantum chemistry log file.
        """
        raise NotImplementedError

    def print_data(self, outfile='NICS.txt'):
        with open(outfile, 'w') as f:
            f.write(self.NICS_data[0].get_header())
            for point in self.NICS_data:
                f.write(point.get_data())

    def ring_analysis(self, ring, tol=0.05, outfile='NICS_ring.txt', lvprt=1):
        """
        Evaluate NICS_iso and NICS_zz for all ghost atoms on the axis of a ring
        and label the NICS(0) and NICS(+-1) points.
        ring    - ring atoms, either as a list of atom indices (starting at 1,
                  in the order of the real atoms in the log file) or as an
                  (n x 3) array of coordinates (Angstrom)
        tol     - maximal distance (Angstrom) of a ghost atom from the ring
                  axis, and from 0 / +-1 Angstrom for the NICS(0) / NICS(1) labels
        outfile - write the table to this file (None: no file)
        Returns a list of dictionaries, one per point on the axis.
        """
        ring = numpy.array(ring)
        if ring.ndim == 1:
            if len(self.atoms) == 0:
                raise ValueError("No atomic coordinates were parsed, give the ring as coordinates")
            for iat in ring:
                if not 1 <= iat <= len(self.atoms):
                    raise ValueError("Ring atom %i not found (%i atoms parsed)"%(iat, len(self.atoms)))
            ring_coords = numpy.array([self.atoms[iat-1]['coor'] for iat in ring])
        else:
            ring_coords = numpy.array(ring, float)

        centroid, normal = ring_frame(ring_coords)

        results = []
        for ipoint, point in enumerate(self.NICS_data):
            rel = point.orig - centroid
            height = numpy.dot(rel, normal)
            offset = numpy.linalg.norm(rel - height * normal)
            if offset > tol:
                continue

            if abs(height) < tol:
                label = 'NICS(0)'
            elif abs(height - 1.) < tol:
                label = 'NICS(1)'
            elif abs(height + 1.) < tol:
                label = 'NICS(-1)'
            else:
                label = 'NICS(%.2f)'%height

            results.append({'point': ipoint, 'height': height, 'offset': offset,
                'NICS_iso': point.ret_NICS_iso(), 'NICS_zz': point.ret_NICS_comp(normal),
                'label': label})

        outstr  = "Ring centroid: (% .5f, % .5f, % .5f)\n"%tuple(centroid)
        outstr += "Ring normal:   (% .5f, % .5f, % .5f)\n"%tuple(normal)
        outstr += "NICS = -sigma (ppm); NICS_zz along the ring normal\n"
        outstr += "%6s %9s %9s %10s %10s  %s\n"%('Point', 'height', 'offset', 'NICS_iso', 'NICS_zz', 'label')
        for res in results:
            outstr += "%6s % 9.4f % 9.4f % 10.4f % 10.4f  %s\n"%('P%i'%res['point'],
                res['height'], res['offset'], res['NICS_iso'], res['NICS_zz'], res['label'])
        if len(results) == 0:
            outstr += " *** No ghost atoms found within %.3f A of the ring axis\n"%tol

        if lvprt >= 1:
            print("\n *** NICS analysis for ring ***")
            print(outstr)
        if outfile is not None:
            with open(outfile, 'w') as f:
                f.write(outstr)
            if lvprt >= 1:
                print("Ring analysis written to %s"%outfile)

        return results

    def vmd_tensors(self, filen='VIST.vmd', vlist=None, scale=1., plot_all=False, lab_min=1000.):
        """
        Draw VMD tensors for all the NICS values parsed.
        filen - name of output file
        vlist - list of dummy atoms to consider
        scale - scale factor for dumb bells
        """
        af = open(filen, 'a')
        self.af = af
        af.write(\
"""axes location Off
display projection Orthographic
display depthcue off
color Display Background white
menu graphics on
mol modstyle 0 %i Licorice 0.100000 30.000000 30.000000
# Uncomment to hide hydrogens
# mol modselect 0 %i not element H

# Use this for rings
mol addrep %i
mol modstyle 1 %i PaperChain 0.050000 10.000000
mol modmaterial 1 %i Glass3
draw delete all

"""%(self.mol_ind, self.mol_ind, self.mol_ind, self.mol_ind, self.mol_ind))

        for ipoint, point in enumerate(self.NICS_data):
            if vlist is None or ipoint in vlist:
                point.vmd_tensor(af, scale=scale, lab_min=lab_min)
                if plot_all:
                    af.write("render TachyonInternal P%i.tga\n"%ipoint)
                    af.write("draw delete all\n")
                af.write("\n")

        af.close()
        print("VMD input appended to %s"%af.name)
        self.mol_ind += 1

class NICS_point:
    """
    Container for NICS data for each individual ghost atom.
    """
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.orig = numpy.array([x, y, z], float)
        self.NICS_iso = None
        self.NICS_tensor = None
        self.evals = None
        self.coor = None

    def __str__(self):
        """
        Print the NICS values using the usual convention that their sign is
        reversed with respect to the shielding tensor.
        """
        outstr = "NICS value at (% .5f, % .5f, % .5f)"%(self.x, self.y, self.z)
        if not self.NICS_iso is None:
            outstr += ": % 6.2f"%(-self.NICS_iso)
        if not self.evals is None:
            outstr += " = ( % 6.2f %+6.2f %+6.2f )/3"%(-self.evals[0], -self.evals[1], -self.evals[2])
        return outstr

    def set_iso(self, NICS_iso):
        """
        Set isotropic NICS value.
        """
        self.NICS_iso = NICS_iso

    def set_tensor(self, NICS_tensor):
        """
        Set NICS tensor.
        """
        self.NICS_tensor = numpy.array(NICS_tensor, float)

    def ret_NICS_iso(self):
        """
        Return the isotropic NICS value (in ppm), i.e. the negative of the
        isotropic shielding.
        """
        if self.NICS_iso is not None:
            return -self.NICS_iso
        return -numpy.trace(self.NICS_tensor) / 3.

    def ret_NICS_comp(self, direction):
        """
        Return the NICS component along the given direction (in ppm):
            NICS_nn = - n^T . sigma . n
        with n the normalised direction. Only the symmetric part of the
        shielding tensor contributes, so the result does not depend on whether
        the program prints sigma or its transpose.
        For the ring normal this is the usual NICS_zz value.
        """
        n = numpy.array(direction, float)
        n /= numpy.linalg.norm(n)
        return -numpy.dot(n, numpy.dot(self.NICS_tensor, n))

    def get_data(self):
        """
        Return data as formatted string.
        """
        outstr = "% 11.7f % 11.7f % 11.7f"%(self.x, self.y, self.z)
        if not self.NICS_iso is None:
            outstr += " % 9.5f"%self.NICS_iso
        if not self.NICS_tensor is None:
            for val in self.NICS_tensor.flatten():
                outstr += " % 9.5f"%val

        return outstr + '\n'

    def get_header(self):
        """
        Return header for output.
        """
        outstr = "%9s  %9s  %9s "%('x', 'y', 'z')
        if not self.NICS_iso is None:
            outstr += "      sigma   "
        if not self.NICS_tensor is None:
            outstr += "   s_XX     s_YX      s_ZX       s_XY      s_YY      s_ZY     s_XZ       s_YZ      s_ZZ"
        return outstr + '\n'

    def diag(self, rightEV=True):
        """
        Diagonalise the NICS tensor.
        """
        # Compute left or right eigenvectors
        if rightEV:
            (evals, coor) = numpy.linalg.eig(self.NICS_tensor)
        else:
            (evals, coor) = numpy.linalg.eig(self.NICS_tensor.T)

        if not numpy.isreal(evals).all():
            with numpy.printoptions(**NICS_PRINTOPTIONS):
                # If there are complex eigenvalues, then these have to be transformed
                #   Form linear combinations of the associated eigenvectors to get a 2-dim irrep
                print("""
*** WARNING: Transforming complex eigenvectors of NICS tensor
      into a real representation""")
                print("NICS Tensor")
                print(self.NICS_tensor)
                print("Eigenvalues")
                print(evals)
                print("Eigenvectors")
                print(coor)

                compl_inds = []
                real_inds  = []
                for i in range(3):
                    if numpy.isreal(evals[i]):
                        real_inds.append(i)
                    else:
                        compl_inds.append(i)
                assert(len(compl_inds)==2)
                assert(numpy.conj(evals[compl_inds[0]]) == evals[compl_inds[1]])

                # Use the real and imaginary parts of the complex-conjugate eigenvectors
                #    to create an invariant subspace
                coor_new = numpy.zeros([3,3])
                coor_new[:, real_inds[0]]  = numpy.real(coor[:, real_inds[0]])
                #coor_new[:, real_inds[0]] /= numpy.linalg.norm(coor_new[:, real_inds[0]]) # already unit "length"

                coor_new[:, compl_inds[0]] = numpy.real(coor[:, compl_inds[0]])
                coor_new[:, compl_inds[0]] /= numpy.linalg.norm(coor_new[:, compl_inds[0]])

                coor_new[:, compl_inds[1]] = numpy.imag(coor[:, compl_inds[0]])
                coor_new[:, compl_inds[1]] /= numpy.linalg.norm(coor_new[:, compl_inds[1]])
                self.coor = coor_new

                print("Real transformation matrix")
                print(self.coor)
                print("Transformed NICS tensor (should be block-diagonal\n with real parts of EVs in diagonal)")
                tmp = numpy.dot(numpy.dot(numpy.linalg.inv(self.coor), self.NICS_tensor), self.coor)
                print(tmp)
                self.evals = numpy.real(evals)
                print("EVs: ", self.evals, "=", numpy.diag(tmp))
        else:
            # eig may return a complex dtype with vanishing imaginary parts
            self.evals = numpy.real(evals)
            self.coor  = numpy.real(coor)

    # VMD part - maybe this should be moved into a separate library
    def vmd_tensor(self, af, scale=1., lab_min=1000.):
        """
        Create a graphical representation of the NICS tensor to be plotted in VMD.
        scale - scale for drawing tensors
        lab_min - minimal value (in ppm) for drawing labels
        """
        self.af = af

        for mu in range(3):
            # The eigenvector is in a column
            #print("A v", mu)
            # Check if the eigenvectors are correct -> yes
              # Av = numpy.dot(self.NICS_tensor, coor[:, mu])
              # print(coor[:, mu], Av, Av / evals[mu])
            #print("Abs")
            #print(abs(coor[:, mu]))

            evmu = self.evals[mu]
            if self.vmd_color(evmu, eps=1.):
                vec = scale * 0.3 * abs(evmu)**.5 * self.coor[:, mu]
                stick = scale * 0.1
                rad = max(scale * 0.03 * abs(evmu)**.5, 0.5 * stick)
                if abs(evmu) >= lab_min:
                    self.plot_quad_comp(self.orig, vec, rad, stick, evmu)
                else:
                    self.plot_quad_comp(self.orig, vec, rad, stick)
                #self.plot_quad_cone(fac, self.orig, coor[:, mu], 0.5 * fac)

    def plot_quad_comp(self, orig, vec, rad, stick=0.1, lab=None):
        """
        Plot component of quadrupole moment.
        """
        self.af.write('draw cylinder ')
        self.vmd_coors(orig - vec)
        self.vmd_coors(orig + vec)
        self.af.write('radius % .3f\n'%stick)

        self.af.write('draw sphere ')
        self.vmd_coors(orig + vec)
        self.af.write('radius % .3f\n'%(2*rad))

        self.af.write('draw sphere ')
        self.vmd_coors(orig - vec)
        self.af.write('radius % .3f\n'%(2*rad))

        if not lab is None:
            self.af.write('draw text ')
            self.vmd_coors(orig + vec * 0.7)
            self.af.write('" % .1f"\n'%lab)

    def plot_quad_cone(self, fac, orig, vec, rad):
        """
        Plot component of quadrupole moment as cone.
        """
        self.af.write('draw cone ')
        self.vmd_coors(orig - fac*vec)
        self.vmd_coors(orig)
        self.af.write('radius % .3f\n'%rad)

        self.af.write('draw cone ')
        self.vmd_coors(orig + fac*vec)
        self.vmd_coors(orig)
        self.af.write('radius % .3f\n'%rad)

    def vmd_coors(self, coor):
        self.af.write('{% .3f % .3f % .3f} '%(coor[0], coor[1], coor[2]))

    def vmd_color(self, val, eps=1.e-3):
        if abs(val) < eps:
            return False
        elif val > 0.0:
            self.af.write('draw color blue\n')
            return True
        else:
            self.af.write('draw color red\n')
            return True

class NICS_parser_gaussian(NICS_parser):
    """
    Parse Gaussian (G09, G16) NICS calculations.

    The ghost atom (Bq) positions are taken from the last orientation table
    printed before the shielding tensors: the "Standard orientation" if
    present, otherwise the "Input orientation". Gaussian reports the shielding
    tensors in the standard orientation unless NoSymm is given, so the
    coordinates and tensors are always in the same frame. (When ghost atoms
    are present, Gaussian usually turns symmetry off and prints only the
    input orientation, so that both frames coincide.)
    Ghost atoms are matched to their tensors through the centre number, which
    also covers Bq atoms given as Z-matrix lines or with freeze codes.
    If the log file contains several shielding sections (e.g. Opt+NMR or
    --Link1-- jobs), only the last one is used.
    If no orientation table is printed, the Cartesian Bq lines of the
    "Symbolic Z-matrix" echo are used as a fallback.
    """
    re_iso  = re.compile(r'^\s*(\d+)\s+(\S+)\s+Isotropic\s*=\s*(\S+)')
    re_comp = re.compile(r'([XYZ]{2})=\s*(-?[0-9.]+(?:[DdEe][-+]?[0-9]+)?)')

    def read(self, logfile, lvprt=1):
        self.NICS_data = []
        self.atoms = []

        last_input = None
        last_standard = None
        zmat_bq = []
        sections = []
        with open(logfile, 'r') as f:
            while True:
                try:
                    line = next(f)
                except StopIteration:
                    break

                if 'Input orientation:' in line:
                    last_input = self.read_orientation(f)
                    last_standard = None
                elif 'Standard orientation:' in line:
                    last_standard = self.read_orientation(f)
                elif 'Symbolic Z-matrix:' in line:
                    zmat_bq = self.read_zmat_bq(f)
                elif 'Magnetic shielding tensor (ppm)' in line:
                    if last_standard is not None:
                        geom, orient = last_standard, 'Standard'
                    elif last_input is not None:
                        geom, orient = last_input, 'Input'
                    else:
                        geom, orient = None, None
                    sections.append({'geom': geom, 'orient': orient,
                        'zmat_bq': list(zmat_bq), 'iso': {}, 'tensor': {}, 'labels': {}})
                elif len(sections) > 0 and 'Isotropic' in line:
                    match = self.re_iso.match(line)
                    if match is None:
                        continue
                    icen, label, iso = int(match.group(1)), match.group(2), float(match.group(3))
                    tensor = []
                    for i in range(3):
                        comps = self.re_comp.findall(next(f))
                        if len(comps) != 3:
                            raise ValueError("Cannot read shielding tensor of centre %i in %s"%(icen, logfile))
                        tensor.append([float(val.replace('D', 'E').replace('d', 'e')) for (lab, val) in comps])
                    sec = sections[-1]
                    sec['iso'][icen] = iso
                    sec['tensor'][icen] = tensor
                    sec['labels'][icen] = label

        if len(sections) == 0:
            raise ValueError("No shielding tensors found in %s"%logfile)
        if len(sections) > 1:
            print(" Note: %i shielding sections found in %s, using the last one"%(len(sections), logfile))
        sec = sections[-1]

        bq_cens = sorted(icen for icen, label in sec['labels'].items() if label == 'Bq')
        if sec['geom'] is not None:
            geom = {row[0]: row for row in sec['geom']}
            if sec['orient'] == 'Standard':
                print(" Using ghost atom coordinates from the Standard orientation of %s"%logfile)
            self.atoms = [{'Z': row[1], 'coor': row[2]} for row in sec['geom'] if row[1] > 0]
            for icen in bq_cens:
                if not icen in geom or geom[icen][1] != 0:
                    raise ValueError("Centre %i of %s is not a ghost atom in the %s orientation"%(icen, logfile, sec['orient']))
                (x, y, z) = geom[icen][2]
                self.NICS_data.append(NICS_point(x, y, z))
        else:
            print(" *** WARNING: no orientation table found in %s,"%logfile)
            print("     using the Bq coordinates of the input echo (correct only with NoSymm)")
            if len(sec['zmat_bq']) != len(bq_cens) or None in sec['zmat_bq']:
                raise ValueError("Cannot determine the Bq coordinates in %s from the input echo"%logfile)
            for (x, y, z) in sec['zmat_bq']:
                self.NICS_data.append(NICS_point(x, y, z))

        for ipoint, icen in enumerate(bq_cens):
            self.NICS_data[ipoint].set_iso(sec['iso'][icen])
            self.NICS_data[ipoint].set_tensor(sec['tensor'][icen])
            self.NICS_data[ipoint].diag()
        print("Finished parsing %s"%logfile)

        if lvprt >= 1:
            print(" *** Printing NICS and eigenvalues ***")
            for ipoint, point in enumerate(self.NICS_data):
                print("P%i ->"%ipoint, point)

    def read_orientation(self, f):
        """
        Read a Gaussian orientation table.
        Returns a list of (centre number, atomic number, numpy.array([x, y, z])).
        """
        for i in range(4): # dashes and the two header lines
            next(f)
        rows = []
        for line in f:
            if '----' in line:
                break
            words = line.split()
            rows.append((int(words[0]), int(words[1]), numpy.array(words[-3:], float)))
        return rows

    def read_zmat_bq(self, f):
        """
        Read the Bq lines of the "Symbolic Z-matrix" echo.
        Cartesian lines with or without a freeze code are converted into
        coordinates, Z-matrix lines are returned as None.
        """
        bq = []
        for line in f:
            words = line.split()
            if len(words) == 0:
                break
            if words[0].split('(')[0] != 'Bq':
                continue
            try:
                if len(words) == 4:
                    bq.append(tuple(float(w) for w in words[1:4]))
                elif len(words) == 5 and words[1] in ('0', '-1'):
                    bq.append(tuple(float(w) for w in words[2:5]))
                else:
                    bq.append(None)
            except ValueError:
                bq.append(None)
        return bq

class NICS_parser_QC(NICS_parser):
    """
    Parse Q-Chem NICS calculations.
    """
    def read(self, logfile, lvprt=2):
        print("Reading Q-Chem information")
        self.NICS_data = []
        with open(logfile, 'r') as f:
            Bqind = 0
            while True:
                try:
                    line = next(f)
                except StopIteration:
                    print("Finished parsing %s"%logfile)
                    break

                if '  GH' in line and not 'atom' in line.lower():
                    words = line.split()
                    if words[1] == 'GH':
                        print(line)
                        (x, y, z) = map(float, words[2:5])
                        self.NICS_data.append(NICS_point(x, y, z))

                if 'ATOM  GH' in line:
                    while True:
                        line = next(f)

                        if 'total shielding tensor' in line:
                            line = next(f)
                            NICS_iso = float(line.split()[-1])
                            self.NICS_data[Bqind].set_iso(NICS_iso)

                            line = next(f)
                            tensor = []
                            for i in range(3):
                                words = next(f).split()
                                tensor.append([float(words[0]), float(words[1]), float(words[2])])
                            self.NICS_data[Bqind].set_tensor(tensor)
                            self.NICS_data[Bqind].diag()

                            Bqind += 1
                            break

        if lvprt >= 1:
            print(" *** Printing NICS and eigenvalues ***")
            for ipoint, point in enumerate(self.NICS_data):
                print("P%i ->"%ipoint, point)

class NICS_parser_TM(NICS_parser):
    """
    Parse TURBOMOLE NICS calculations.
    """
    def read(self, logfile, lvprt=2):
        print("Reading TURBOMOLE information")
        self.NICS_data = []
        with open(logfile, 'r') as f:
            Bqind = 0
            while True:
                try:
                    line = next(f)
                except StopIteration:
                    print("Finished parsing %s"%logfile)
                    break

                if '  atomic coordinates  ' in line:
                    while True:
                        line = next(f)
                        words = line.split()
                        if not len(words) == 6:
                            break
                        el = words[3]
                        if el == 'q':
                            (x, y, z) = map(float, line.split()[0:3])
                            x *= units.length['A']
                            y *= units.length['A']
                            z *= units.length['A']
                            self.NICS_data.append(NICS_point(x, y, z))

                if 'ATOM  q' in line:
                    while True:
                        line = next(f)

                        if 'total magnetic shielding' in line:
                            line = next(f)
                            NICS_iso = float(line.split()[-1])
                            self.NICS_data[Bqind].set_iso(NICS_iso)

                            line = next(f)
                            tensor = []
                            for i in range(3):
                                words = next(f).split()
                                tensor.append([float(words[0]), float(words[1]), float(words[2])])
                            self.NICS_data[Bqind].set_tensor(tensor)
                            self.NICS_data[Bqind].diag()

                            Bqind += 1
                            break

        if lvprt >= 1:
            print(" *** Printing NICS and eigenvalues ***")
            for ipoint, point in enumerate(self.NICS_data):
                print("P%i ->"%ipoint, point)

class NICS_parser_dalton(NICS_parser):
    """
    Parse Dalton NICS calculations.
    """
    def read(self, logfile, lvprt=1):
        print("Reading Dalton information")
        self.NICS_data = []
        with open(logfile, 'r') as f:
            bq_coords = []
            for line in f:
                if 'Cartesian Coordinates (xyz format; angstrom)' in line:
                    in_coords = False

                    while True:
                        try:
                            ln = next(f)
                        except StopIteration:
                            break

                        words = ln.split()

                        if not words and not in_coords:
                            continue
                        if not words and in_coords:
                            break
                        if len(words) < 4:
                            continue
                        try:
                            x = float(words[1])
                            y = float(words[2])
                            z = float(words[3])
                        except ValueError:
                            continue

                        in_coords = True

                        if words[0] == 'Bq':
                            bq_coords.append((x, y, z))

        if not bq_coords:
            raise RuntimeError(
                "No Bq coordinates found in the 'Cartesian Coordinates "
                "(xyz format; angstrom)' section of %s" % logfile
            )
        for (x, y, z) in bq_coords:
            self.NICS_data.append(NICS_point(x, y, z))

        with open(logfile, 'r') as f:
            Bqind = 0
            while True:
                try:
                    line = next(f)
                except StopIteration:
                    print("Finished parsing %s" % logfile)
                    break

                if 'Chemical shielding for Bq' in line:
                    while True:
                        line = next(f)
                        if 'Shielding constant' in line:
                            break

                    NICS_iso = float(line.split()[2])

                    while True:
                        line = next(f)
                        if 'Total shielding tensor' in line:
                            break

                    tensor = []

                    while len(tensor) < 3:
                        line = next(f).strip()
                        if not line:
                            continue

                        words = line.split()

                        if not words or words[0] != 'Bq':
                            continue

                        #print("DEBUG:", words)

                        try:
                            vals = [float(x) for x in words[-3:]]
                        except ValueError:
                            raise RuntimeError(f"Malformed tensor row: {words}")

                        tensor.append(vals)

                        #tensor.append([
                            #float(words[3]),
                            #float(words[4]),
                            #float(words[5])
                        #])

                    if Bqind >= len(self.NICS_data):
                        raise RuntimeError(
                            "More shielding tensors than Bq coordinates found"
                        )

                    self.NICS_data[Bqind].set_iso(NICS_iso)
                    self.NICS_data[Bqind].set_tensor(tensor)
                    self.NICS_data[Bqind].diag()
                    Bqind += 1

        if lvprt >= 1:
            print(" *** Printing NICS and eigenvalues ***")
            for ipoint, point in enumerate(self.NICS_data):
                print("P%i ->" % ipoint, point)
