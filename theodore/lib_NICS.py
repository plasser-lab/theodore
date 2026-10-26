"""
Analysis routines in the context of NICS (nucleus independent chemical shift)
calculations.

Author: Felix Plasser
Date: Nov. 2020
"""

import numpy
from theodore import units

# Input


# Output

class NICS_parser:
    """
    Parse and visualise NICS calculations.
    """
    def __init__(self, logfile=None):
        self.NICS_data = []
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
            self.evals = evals
            self.coor  = coor

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

class NICS_parser_g09(NICS_parser):
    """
    Parse Gaussian NICS calculations.
    """
    def read(self, logfile, lvprt=1):
        self.NICS_data = []
        with open(logfile, 'r') as f:
            Bqind = 0
            while True:
                try:
                    line = next(f)
                except StopIteration:
                    print("Finished parsing %s"%logfile)
                    break

                if line[0:4] == ' Bq ':
                    (x, y, z) = map(float, line.split()[1:4])
                    self.NICS_data.append(NICS_point(x, y, z))

                if 'Bq   Isotropic' in line:
                    NICS_iso = float(line.split()[4])
                    self.NICS_data[Bqind].set_iso(NICS_iso)

                    tensor = []
                    for i in range(3):
                        words = next(f).split()
                        tensor.append([float(words[1]), float(words[3]), float(words[5])])
                    self.NICS_data[Bqind].set_tensor(tensor)
                    self.NICS_data[Bqind].diag()

                    Bqind += 1

        if lvprt >= 1:
            print(" *** Printing NICS and eigenvalues ***")
            for ipoint, point in enumerate(self.NICS_data):
                print("P%i ->"%ipoint, point)

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
