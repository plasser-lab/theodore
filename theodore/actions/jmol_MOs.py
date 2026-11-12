#!/usr/bin/env python3
"""
Automatic plotting of MOs with jmol.
"""

from __future__ import print_function, division
from .actions import Action
import sys, os
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    lib_mo = importer.lazy_import_as('..lib_mo', 'lib_mo')
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')
    lib_file = importer.lazy_import_as('..lib_file', 'lib_file')
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    input_options = importer.lazy_import_as('..input_options', 'input_options')


class mo_output:
    """
    Abstract base class for MO output.
    """
    def __init__(self, moc, jopt):
        self.moc = moc
        self.jopt = jopt
        self.outstr = ''

    def output(self, ofileh):
        self.pre()
        self.print_mos()
        self.post(ofileh)

    def mopath(self, mo, of='png'):
        if of == 'pngt': of = 'png'
        return "MO_%s.%s"%(mo,of)

    def pre(self):
        raise error_handler.PureVirtualError()

    def print_mos(self):
        raise error_handler.PureVirtualError()

    def post(self, ofileh):
        ofileh.write(self.outstr)

class mo_output_jmol(mo_output):
    """
    MO output in standard jmol format.
    """
    def pre(self):
        if self.moc.mldfile != "":
            self.outstr += '\nload "' + self.moc.mldfile + '" FILTER "nosort"\n'
        self.outstr += "mo titleformat ''\n"
        if self.jopt['rot_best']:
            self.outstr += "rotate best\n"
        if self.jopt['rot_custom']:
            self.outstr += "rotate x %.3f\n"%self.jopt['rot_x']
            self.outstr += "rotate y %.3f\n"%self.jopt['rot_y']
            self.outstr += "rotate z %.3f\n"%self.jopt['rot_z']
        self.outstr += "background white\n" + "mo fill\n"
        self.outstr += "mo cutoff %.3f\n\n"%self.jopt['cutoff']

    def print_mos(self):
        squ_string_plus = ''
        squ_string_minus = ''
        for imo in self.moc.molist:
            if self.jopt['ccode']:
                if self.moc.is_occ(imo):
                    self.outstr += "mo color %s\n"%self.jopt['chole']
                else:
                    self.outstr += "mo color %s\n"%self.jopt['celec']
            self.outstr += "mo %s\n"%self.moc.moname(imo)
            self.outstr += "write image %s \"%s\"\n"%(self.jopt['oformat'], self.moc.mopath(imo, self.jopt['oformat']))

            if self.jopt['do_dens']:
                occ = self.moc.ret_occ(imo)
                # One has to take the squareroot here for jmol to undlerstand this correctly
                if occ > 0.:
                    squ_string_plus += ' %.3f %i'%(occ**.5, imo+1)
                if occ < 0.:
                    squ_string_minus += ' %.3f %i'%((-occ)**.5, imo+1)

        if squ_string_plus != '':
            self.outstr += "mo color %s\n"%self.jopt['celec']
            self.outstr += 'mo [%s] squared\n'%squ_string_plus
            self.outstr += "write image %s \"%s\"\n"%(self.jopt['oformat'], self.moc.denspath("plus", self.jopt['oformat']))
        if squ_string_minus != '':
            self.outstr += "mo color %s\n"%self.jopt['chole']
            self.outstr += 'mo [%s] squared\n'%squ_string_minus
            self.outstr += "write image %s \"%s\"\n"%(self.jopt['oformat'], self.moc.denspath("minus", self.jopt['oformat']))

class mo_output_html(mo_output):
    """
    HTML file for visualizing the MOs created with jmol.
    """
    def pre(self):
        self.htable = lib_file.htmltable(ncol=2)

    def print_mos(self):
        for imo in self.moc.molist:
            el = '<img src="%s" "border="1" width="%i">'%(self.moc.mopath(imo, self.jopt['oformat']), self.jopt['width'])
            el += self.moc.mo_extra(imo,pref="<br> MO %i:"%(imo+1))
            self.htable.add_el(el)

        el = '<img src="%s" "border="1" width="%i">'%(self.moc.denspath("plus", self.jopt['oformat']), self.jopt['width'])
        el += "<br> Density (plus)"
        self.htable.add_el(el)
        el = '<img src="%s" "border="1" width="%i">'%(self.moc.denspath("minus", self.jopt['oformat']), self.jopt['width'])
        el += "<br> Density (minus)"
        self.htable.add_el(el)

    def post(self, ofileh):
        ofileh.write("<h2>Orbitals - %s</h2>\n"%self.moc.mldfile)
        self.htable.close_table()
        ofileh.write(self.htable.ret_table())

class mo_output_tex(mo_output):
    """
    tex file for visualizing the MOs created with jmol.
    """
    def pre(self):
        self.ltable = lib_file.latextabular(ncol=2)

    def print_mos(self):
        moex = []
        for imo in self.moc.molist:
            el = "\\incMO{%s}"%(self.moc.mopath(imo, self.jopt['oformat']))
            moex += [self.moc.mo_extra(imo)]
            lastcol = self.ltable.add_el(el)
            if lastcol:
                self.ltable.add_row(moex)
                moex = []

    def post(self, ofileh):
        #self.ltable.close_table() # not needed??
        ofileh.write(self.ltable.ret_table())

class mocoll:
    """
    MO output for start and end indices.
    Also base class for other options.
    """
    def __init__(self, st_ind, en_ind, mldfile=""):
        self.mldfile = mldfile

        if not mldfile=="":
            self.moset = lib_mo.MO_set_molden(mldfile)
            self.moset.read(lvprt=0)
            maxmo = self.moset.ret_num_mo()
        else:
            maxmo = 10000

        self.molist = list(range(st_ind-1, min(en_ind, maxmo)))
        self.virtlist = []

    def moname(self, imo):
        return str(imo+1)

    def is_occ(self, imo):
        """
        Return True if the MO is occupied, and False if not.
        """
        return True

    def mopath(self, imo, of='png'):
        if of == 'pngt': of = 'png'
        return "%sMO_%s.%s"%(self.ret_label(),self.moname(imo),of)

    def denspath(self, dname, of='png'):
        if of == 'pngt': of = 'png'
        return "%sdens_%s.%s"%(self.ret_label(),dname,of)

    def mo_extra(self, imo, pref="", postf=""):
        if self.mldfile=="":
            return ""
        else:
            sym = self.moset.ret_sym(imo)
            try:
                ene, occ = self.moset.ret_eo(imo)
            except:
                print("\n ERROR: imo = %i"%imo)
                print("ens:", self.moset.ens)
                print("occs:", self.moset.occs)
                raise
            return "%s %5s %.4f / %.4f %s"%(pref,sym,ene,occ,postf)

    def ret_occ(self, imo):
        """
        Return occupancy of orbital (including eneocc option).
        """
        ene, occ = self.moset.ret_eo(imo)
        if self.eneocc:
            return ene
        else:
            return occ

    def ret_label(self):
        if self.mldfile=="":
            return ""
        else:
            #return self.mldfile.split('.')[0] + '_'
            return self.mldfile.replace('.','-') + '_'

class mocollf(mocoll):
    """
    For frontier MOs.
    """
    def __init__(self, en_ind, mldfile=""):
        mocoll.__init__(self, 1, 2*en_ind, mldfile)

    def moname(self, imo):
        if imo==0:
            return "homo"
        elif imo==1:
            return "lumo"
        elif imo%2==0:
            return "homo-%i"%(imo/2)
        else:
            return "lumo+%i"%(imo/2)

    def is_occ(self, imo):
        return imo%2==0

    def mo_extra(self, imo, pref="", postf=""):
        if self.mldfile=="":
            return ""
        else:
            ihomo = self.moset.ret_ihomo()
            if imo%2==0:
                imo2 = ihomo - imo//2
            else:
                imo2 = ihomo + imo//2 + 1
            try:
                ene, occ = self.moset.ret_eo(imo2)
            except:
                print("\n ERROR: imo2 = %i"%imo2)
                raise
            return "%s %.3f / %.3f %s"%(pref,ene,occ,postf)

class mocoll_occ(mocoll):
    """
    Specify occupation threshold.
    """
    def __init__(self, occmin=0.01, occmax=2.0, mldfile="", eneocc=False):
        self.eneocc = eneocc
        if mldfile == "":
            raise error_handler.MsgError("mldfile has to be specified for occupation screening!")

        self.mldfile = mldfile
        self.moset = lib_mo.MO_set_molden(mldfile)
        self.moset.read(lvprt=0)

        self.molist = []
        if eneocc:
            loccs = self.moset.ens
        else:
            loccs = self.moset.occs
        for imo, occ in enumerate(loccs):
            if occmin <= abs(occ) <= occmax:
                self.molist.append(imo)

    def is_occ(self, imo):
        """
        Check if orbital is occupied. This only works if the
        occupied (hole) orbitals are given with negative populations.
        """
        if self.eneocc:
            occnum = self.moset.ens[imo]
        else:
            occnum = self.moset.occs[imo]
        return occnum < 0 or occnum > 1

class jmol_options(input_options.write_options):
    def jmol_input(self, nfiles=1):
        # some defaults
        self['do_dens'] = False
        self['chole']   = 'blue red'
        self['celec']   = 'red blue'

        self.read_float('Cutoff value', 'cutoff', 0.05)

        #print ""
        self.choose_list(
            'Specification of the orbital indices to be plotted',
            'spec',
        [
            ('sten', 'Start and end indices'),
            ('frontier', 'Number of frontier orbitals'),
            ('occ', 'Occupation threshold')
        ]
        )
        #self.read_yn('Specification in terms of frontier orbitals', 'fr_mos')

        if self['spec'] == 'sten':
            self.read_int('First orbital index to be plotted', 'st_ind', 1)
            self.read_int('Last orbital index to be plotted',  'en_ind', 10)
            self['preprocess'] = False
        elif self['spec'] == 'frontier':
            self.read_int('Number of frontier orbitals',  'en_ind', 3)
            self['preprocess'] = False
        elif self['spec'] == 'occ':
            self.read_float('Minimal absolute occupancy', 'occmin', 0.01)
            self.read_float('Maximal absolute occupancy', 'occmax', 1.99)

            self.read_yn('Preprocess and merge the Molden files', 'preprocess', nfiles>1)
            if not self['preprocess']:
                self.read_yn('Compute occupancy-weighted densities (attachment/dentachment, hole/electron, etc)?', 'do_dens', False)

            self.read_yn('Interpret energies as occupations (use for Q-Chem)', 'eneocc', False)

        else:
            raise error_handler.ElseError(self['spec'], 'spec')

        self.read_yn('Use "rotate best" command (available since Jmol 14)', 'rot_best', True)
        self.read_yn('Additional custom rotation of the molecule?', 'rot_custom')
        if self['rot_custom']:
            self.read_float('Rotation around the x-axis', 'rot_x', 0.)
            self.read_float('Rotation around the y-axis', 'rot_y', 0.)
            self.read_float('Rotation around the z-axis', 'rot_z', 0.)

        self.read_yn('Color code for hole / electron orbitals', 'ccode', True)
        if self['ccode']:
            self.read_str('Hole (occupied) orbitals', 'chole', 'blue red')
            self.read_str('Electron (virtual) orbitals', 'celec', 'orange green')
        self.read_str('Format of the output files (png, pngt, ...)', 'oformat', 'pngt')
        self.read_int('Width of images in output html file', 'width', 400)

        self.read_yn('Run Jmol?', 'run_jmol', False)

    def preprocess(self, mldfiles, out='merged.mld'):
        f = open(out, 'w')
        f.write('[Molden Format]\n')

        mos = lib_mo.MO_set_molden(mldfiles[0])
        mos.read()
        for line in mos.header.split('\n')[1:]:
            f.write(line+'\n')
        f.write('[MO]\n')
        f.write(mos.ret_coeffs(self['occmin'], self['occmax'], self['eneocc'], sym=mldfiles[0]))

        for mldfile in mldfiles[1:]:
            mos = lib_mo.MO_set_molden(mldfile)
            mos.read()
            f.write(mos.ret_coeffs(self['occmin'], self['occmax'], self['eneocc'], sym=mldfile))

        f.close()



class JMolMOs(Action):

    name = 'jmol_mos'

    _colt_description = 'Orbital/density plotting in Jmol'

    _user_input = """
    # List of Molden files with orbitals
    mldfiles = :: list(str), optional, alias=f
    """

    _lazy_imports = LazyImporter({
            '..lib_mo': 'lib_mo',
            '..theo_header': 'theo_header',
            '..error_handler': 'error_handler',
            '..input_options': 'input_options',
            '..lib_file': 'lib_file',
    })

    def run(mldfiles):
        theo_header.print_header(__class__._colt_description)

        if mldfiles is None:
            print("No file specified, generating generic script")
            mldfiles = ['']
            pref = ''
        elif len(mldfiles) == 1:
            print("Analyzing the file:", mldfiles[0])
            pref = mldfiles[0] + '.'
        else:
            print("Analyzing the files:", mldfiles)
            pref = 'multi.'

        jopt = jmol_options('jmol.in')
        jopt.jmol_input()

        if jopt['preprocess']:
            jopt.preprocess(mldfiles)
            mldfiles = ['merged.mld']

        jo = lib_file.wfile('jmol_orbitals.spt')
        ho = lib_file.htmlfile('%sorbitals.html'%pref)
        lo = lib_file.latexfile('%sorbitals.tex'%pref)

        ho.pre('Orbitals')
        lo.pre(None, graphicx=True, docclass='{standalone}')

        for mldfile in mldfiles:
            print('Analyzing %s ...\n'%mldfile)
            if jopt['spec'] == 'sten':
                moc = mocoll(jopt['st_ind'], jopt['en_ind'], mldfile)
            elif jopt['spec'] == 'frontier':
                moc = mocollf(jopt['en_ind'], mldfile)
            elif jopt['spec'] == 'occ':
                moc = mocoll_occ(jopt['occmin'], jopt['occmax'], mldfile, jopt['eneocc'])
            else:
                raise error_handler.ElseError(self['spec'], 'spec')

            moout = mo_output_jmol(moc, jopt)
            moout.output(jo)

            moh = mo_output_html(moc, jopt)
            moh.output(ho)

            mol = mo_output_tex(moc, jopt)
            mol.output(lo)

        ho.post(lvprt=1)
        print("  -> View in browser.")
        lo.post(lvprt=1)
        print("  -> Compile with pdflatex (or adjust first).")

        jo.post(lvprt=1)
        if mldfiles == [""]:
            print("  -> Open the Molden-file in jmol and execute the commands contained in this file.")
        else:
            if jopt['run_jmol']:
                import subprocess
                print("Running jmol ...")

                subprocess.call(["jmol", "-n", jo.name])
            else:
                print("  -> Now simply run \"jmol -n %s\" to plot all the orbitals.\n"%jo)
