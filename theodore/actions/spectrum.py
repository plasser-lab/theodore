#!/usr/bin/env python3
"""
Create a convoluted spectrum from the oscillator strengths.
"""

from __future__ import print_function, division
import math, numpy
from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    units = importer.lazy_import_as('..units', 'units')
    lib_file = importer.lazy_import_as('..lib_file', 'lib_file')
    input_options = importer.lazy_import_as('..input_options', 'input_options')
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')
    matplotlib = importer.lazy_import('matplotlib')
    pylab = importer.lazy_import('pylab')


do_plots = True

class spec_options(input_options.write_options):
    def spec_input(self):
        self.read_int("Number of points in the spectrum", "npts", 200)
        self.read_float("Minimum energy in the plot (eV)", "emin", 2.0)
        self.read_float("Maximum energy in the plot (eV)", "emax", 8.0)
        self.read_float("FWHM broadening (eV)", "fwhm", 0.5)

        self.choose_list(
            "What do you want to plot?",
            "weight",
        [
            (1, "Spectrum weighted by oscillator strength"),
            (2, "Density of states (uniform weighting)"),
            (3, "Absorption cross section"),
            (4, "Molar extinction coefficient"),
            (5, "Radiative emission rate")
        ], 1)

        self.read_int("Lineshape: 1 - Lorentzian, 2 - Gaussian", "lineshape", 1)
        self.spec = spectrum(**self.opt_dict)

        self.read_yn("Use restrictions?", "restr", False)
        if self['restr']:
            self.restrictions()

        if self["weight"] <= 2:
            self.read_yn("Normalize the spectrum?", "normalize", False)
        else:
            self["normalize"] = False

    def restrictions(self):
        self['rlist'] = []

        for irestr in range(1, 1000):
            rstr = self.ret_str('Input restriction #%i (e.g. "CT > 0.5")'%irestr)
            if rstr == '': break

            words = rstr.split()
            self['rlist'].append((words[0], words[1], words[2]))

        print()

    def make_spec(self, lvprt=2, ylabel=None):
        for filen in self['ana_files']:
            sfile = lib_file.summ_file(filen)
            header = sfile.ret_header()
            ddict  = sfile.ret_ddict()
            state_labels = sfile.ret_state_labels()

            n = ntake = 0

            for state in state_labels:
                n += 1
                try:
                    f = ddict[state]['f']
                except KeyError:
                    f = 0.

                if self['restr']:
                    for restr in self['rlist']:
                        evalstr = ("ddict[state]['%s'] %s %s")%restr
                        if not eval(evalstr):
                            break
                    else:
                        self.spec.add(f, ddict[state]['dE(eV)'])
                        ntake += 1
                else:
                    self.spec.add(f, ddict[state]['dE(eV)'])
                    ntake += 1

            if lvprt >= 2:
                print("Considering %3i out of %3i states from %s"%(ntake, n, filen))

        if lvprt >= 1:
            self.spec.info()

        #if self['normalize']:
            #self.spec.normalize()

        self.spec.ascii_file()

        if do_plots:
            self.spec.plot(xunit='eV',  weight=self['weight'], normalize=self['normalize'], ylabel=ylabel)
            self.spec.plot(xunit='nm',  weight=self['weight'], normalize=self['normalize'], ylabel=ylabel)
            self.spec.plot(xunit='rcm', weight=self['weight'], normalize=self['normalize'], ylabel=ylabel)

class gauss:
    def __init__(self,fwhm):
        """
        Initialise spectrum; fwhm in eV.
        """
        self.sig = fwhm * (8 * math.log(2))**(-0.5)
        self.N = units.energy['eV'] / self.sig * (2*numpy.pi)**(-.5)

    def ev(self,A,x0,x):
        """
        Evaluate spectrum at point x.
        """
        return self.N * A * math.exp( -0.5*((x-x0)/self.sig)**2)

class lorentz:
    def __init__(self,fwhm):
        """
        Initialise spectrum; fwhm in eV.
        """
        self.c = 0.25*fwhm**2
        self.N = 2 / numpy.pi * units.energy['eV'] / fwhm # Normalisation factor

    def ev(self,A,x0,x):
        """
        Evaluate spectrum at point x.
        """
        return self.N * A/( (x-x0)**2/self.c + 1)

class spectrum:
    def __init__(self,npts,emin,emax,fwhm,weight,lineshape,ana_files):
      self.npts=npts
      (self.emin, self.emax) = (emin, emax)
      if lineshape==1:
          self.f=lorentz(fwhm)
      elif lineshape==2:
          self.f=gauss(fwhm)
      self.weight = weight

      self.en = [emin + float(i)/self.npts*(emax-emin) for i in range(self.npts+1)]       # the energy grid needs to be calculated only once
      lamev = units.energy['nm'] * units.energy['eV']
      self.lam = [lamev / en for en in self.en]
      self.spec=numpy.array([ 0. for i in range(self.npts+1) ])
      self.dos=numpy.array([ 0. for i in range(self.npts+1) ])

      self.sticks=[] # list of pairs (A,x0), unit for x0: eV

    def add(self,A,x0):
        for i in range(self.npts+1):
            self.dos[i] +=self.f.ev(.1,x0,self.en[i])

        if A!=0.:
            pre = 1
            if self.weight == 1: # Oscillator strength
                pre = 1 / self.f.N # Take out normalisation factor to make sticks as big as the oscillator strength
            elif self.weight == 3: # Absorption cross section
                pre = 2 * numpy.pi**2 / units.constants['c0'] * units.length['A']**2
                # Einstein B. Does not look correct
                #B   = 2 * numpy.pi**2 / units.constants['c0'] / (x0 / units.energy['eV']) * units.time['s'] / units.mass['kg']
                #print('B (Einstein) at %.6f: %.6f s/kg'%(x0, B*A))
            elif self.weight == 4: # Exctinction coefficient
                pre  = 2 * numpy.pi**2 / units.constants['c0'] * units.length['cm']**2
                pre *= units.constants['Nl'] / numpy.log(10) / 1000
                pre *= 1e-5
            elif self.weight == 5: # Differential emission rate
                pre  = 2 / units.constants['c0']**3
                pre *= (x0 / units.energy['eV'])**2
                AA   = pre / units.time['s']
                print('A (Einstein) at %.4f: %12.2f s-1'%(x0, AA*A))
            self.sticks += [(self.f.N * pre * A,x0)]
            for i in range(self.npts+1):
                self.spec[i]+=self.f.ev(pre*A,x0,self.en[i])

    def info(self):
        print("\nSpectrum constructed from %i states with non-vanishing osc. strength"%len(self.sticks))

    def normalize(self):
        """
        Note: this is deactivated.
        """
        smax = max(self.spec)
        print('Normalizing the spectrum...')
        print('Maximum: % .5f'%smax)
        for i, sp in enumerate(self.spec):
            self.spec[i] = sp / smax

        smax = max(self.dos)
        print('Normalizing the DOS ...')
        print('Maximum: % .5f'%smax)
        for i, sp in enumerate(self.dos):
            self.dos[i] = sp / smax

        #for i, stick in enumerate(self.sticks):
         #   self.sticks[i] = (stick[0] / smax, stick[1])

    def ascii_file(self, fname='spectrum.dat'):
        wf = lib_file.wfile(fname)
        wf.write('   eV     spectrum    DOS     nm\n')

        wt = lib_file.asciitable(ncol=4)
        for i, en in enumerate(self.en):
            wt.add_row([en, self.spec[i], self.dos[i], self.lam[i]])

        wf.write(wt.ret_table())
        wf.post(lvprt=1)

    def plot(self, xunit=1, lvprt=1, weight=1, normalize=True, ylabel=None):
        pylab.figure(figsize=(10,7))
        matplotlib.rc('font', size=22)

        if xunit.lower() == 'ev':
            xlist = self.en
            pylab.xlabel('Energy (eV)')
            plot_sticks = self.sticks
            (xmin, xmax) = (self.emin, self.emax)
            pname = 'ev.png'
        elif xunit.lower() == 'nm':
            xfac = units.energy['nm'] * units.energy['eV']
            xlist = self.lam
            #pylab.xlabel(r'$\lambda$') not working ...
            pylab.xlabel('Wavelength (nm)')
            plot_sticks = [(A, xfac/x0) for A,x0 in self.sticks]
            (xmin, xmax) = (xfac / self.emax, xfac / self.emin) # reverse the plot
            pname = 'nm.png'
        elif xunit.lower() == 'rcm':
            xfac = 1. / units.energy['eV'] * units.energy['rcm']
            xlist = [en * xfac for en in self.en]
            pylab.xlabel('Wavenumber (1/cm)')
            plot_sticks = [(A, x0 * xfac) for A,x0 in self.sticks]
            (xmin, xmax) = (self.emin * xfac, self.emax * xfac)
            pname = 'rcm.png'
        else:
            raise error_handler.ElseError('xunit', xunit)

        if weight == 2: # DOS
            if ylabel is None:
                ylabel = "DOS"
            if normalize:
                pylab.plot(xlist, self.dos / max(self.dos), 'k-')
                ymax = 1
            else:
                pylab.plot(xlist, self.dos, 'k-')
                ymax = max(self.dos)
            pname = 'dos_' + pname
        else:
            if weight == 1: # oscillator strength
                if ylabel is None:
                    ylabel = "Oscillator strength"
                pname = 'f_' + pname
            elif weight == 3: # Absorption cross section
                ylabel = "Absorption cross section (Ang^2)"
                pname = 'cross_' + pname
            elif weight == 4: # Extinction coefficient
                ylabel = r'$\epsilon/10^5$ (M$^{-1}$cm$^{-1}$)'
                pname = 'eps_' + pname
            elif weight == 5: # Emission rate
                ylabel = "Emission rate (dim. less)"
                pname = 'emission_' + pname
            else:
                raise error_handler.ElseError('weight', weight)
            if normalize:
                pylab.plot(xlist, self.spec / max(self.spec), 'k-')
                ymax = 1.
            else:
                pylab.plot(xlist, self.spec, 'k-')
                ymax = max(self.spec)
            for A,x0 in plot_sticks:
                pylab.plot([x0, x0], [-1., A], 'rx-')

        pylab.ylabel(ylabel)
        pylab.axis(xmin=xmin, xmax=xmax, ymin=0., ymax=ymax)
        pylab.savefig(pname, dpi=300)

        if lvprt >= 1:
            print("Spectrum file %s created."%pname)

class Spectrum(Action):

    name = 'spectrum'

    _colt_description = 'Convoluted spectrum from analyze_tden output'

    _user_input = """
    # Files produced by analyze_tden.py
    tden_summs = :: list(existing_file)
    """

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..units': 'units',
            '..lib_file': 'lib_file',
            '..input_options': 'input_options',
            '..error_handler': 'error_handler',
            'matplotlib': 'matplotlib',
            'pylab': 'pylab',
    })

    def run(tden_summs):
        matplotlib.use('Agg')
        theo_header.print_header(title=__class__._colt_description)

        sopt = spec_options('spectrum.in')
        sopt['ana_files'] = tden_summs
        sopt.spec_input()

        sopt.make_spec()
