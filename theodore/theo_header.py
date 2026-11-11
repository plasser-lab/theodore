from __future__ import print_function, division
import os

width=80

def print_header(*args, **kwargs):
    print((ret_header(*args, **kwargs)))

def ret_header(title=None, ioptions=None, cfile=None, ver='3.2.1'):
    hstr  = width*'=' + '\n'

    hstr += addlinec("TheoDORE %s"%ver)
    hstr += addlinec("Theoretical Density, Orbital Relaxation and Exciton analysis")
    hstr += addlinec()
    hstr += addlinec("Author: Felix Plasser")
    hstr += addlinec("Contributions by: L. Stojanovic, G. Hermann, S. Mai,")
    hstr += addlinec(" M.F.S.J. Menger, P. Kimber, D. Xing, S. Ghosh, R. G. Bahumanya")

    hstr += width*'-' + '\n'

    hstr += addlinec("References for the modules used")
    hstr += addlinec("(see also http://theodore-qc.sourceforge.net/literature.html)")

    hstr += add_stden(cfile)
    hstr += add_exciton(ioptions)
    hstr += add_extra(ioptions)
    hstr += add_cclib(ioptions)
    hstr += add_orbkit(ioptions)
    hstr += add_VIST(cfile)

    hstr += addlinec()
    hstr += addlinel("Program citation:", 3)
    hstr += addlinel("F. Plasser, J. Chem. Phys. (2020), 152, 084108.")

    if not title==None:
        hstr += width*'-' + '\n'
        hstr += addlinec(title)

    hstr += width*'=' + '\n'

    return hstr

def addlinec(line=""):
    return "|" + line.center(width-2) + "|\n"

def addlinel(line="", lpad=5):
    return "|" + lpad*' ' + line.ljust(width-2-lpad) + "|\n"

def add_stden(cfile):
    if not cfile is None:
        cfileb = os.path.basename(cfile)
    else:
        return ''

    rstr = ''
    if 'analyze' in cfileb:
        rstr += addlinec()
        rstr += addlinel("Transition density matrix analysis:", 3)
        rstr += addlinel("F. Plasser and H. Lischka")
        rstr += addlinel("J. Chem. Theory Comput. (2012), 8, 2777.")
        rstr += addlinec()
        rstr += addlinel("Transition and difference density matrix analysis:", 3)
        rstr += addlinel("F. Plasser, M. Wormit, A. Dreuw")
        rstr += addlinel("J. Chem. Phys. (2014), 141, 024106.")

    return rstr

def add_exciton(ioptions):
    try:
        prop_list = ioptions['prop_list']
    except TypeError:
        return ''

    rstr = ''

    if ('RMSeh' in prop_list) or ('dexc' in prop_list):
        rstr += addlinec()
        rstr += addlinel("Exciton analysis:", 3)
        rstr += addlinel("S. A. Baeppler, F. Plasser, M. Wormit, A. Dreuw")
        rstr += addlinel("Phys. Rev. A (2014), 90, 052521.")

    if 'RMSeh' in prop_list:
        rstr += addlinec()
        rstr += addlinel("Approximate RMSeh/dexc formula:", 3)
        rstr += addlinel("S. A. Mewes, J.-M. Mewes, A. Dreuw, F. Plasser")
        rstr += addlinel("Phys. Chem. Chem. Phys. (2016), 18, 2548.")

    if 'dH-E' in prop_list or 'Corr' in prop_list or 'sigH' in prop_list or 'sigE' in prop_list:
        rstr += addlinec()
        rstr += addlinel("Statistical analysis of excitations:", 3)
        rstr += addlinel("F. Plasser, B. Thomitzni, S. A. Baeppler et al.")
        rstr += addlinel("J. Comput. Chem. (2015), 36, 1609.")

    if ioptions['comp_dntos']:
        rstr += addlinec()
        rstr += addlinel("Conditional densities and DNTOs:", 3)
        rstr += addlinel("F. Plasser")
        rstr += addlinel("ChemPhotoChem (2019), 3, 702.")

    return rstr

def add_extra(ioptions):
    try:
        prop_list = ioptions['prop_list']
    except TypeError:
        return ''

    rstr = ''

    if ('S_HE' in prop_list) or ('Z_HE' in prop_list):
        rstr += addlinec()
        rstr += addlinel("Electron-hole entanglement:", 3)
        rstr += addlinel("F. Plasser")
        rstr += addlinel("J. Chem. Phys. (2016), 144, 194107.")

    if ('QTa' in prop_list) or ('LOC' in prop_list):
        rstr += addlinec()
        rstr += addlinel("Analysis of ionic states:", 3)
        rstr += addlinel("S. A. do Monte, R. F. K. Spada, et al.")
        rstr += addlinel("J. Phys. Chem. A, (2023), 127, 9842.")

    return rstr

def add_cclib(ioptions):
    try:
        rtype = ioptions['rtype'].lower()
    except TypeError:
        return ''

    rstr = ''

    if rtype in ['cclib', 'gamess', 'orca']:
        rstr += addlinec()
        rstr += addlinel("cclib for structure parsing (http://cclib.github.io):", 3)
        rstr += addlinel("N. M. O'Boyle, A. L. Tenderholt, K. M. Langner")
        rstr += addlinel("J. Comput. Chem. (2008), 29, 839.")

    return rstr

def add_orbkit(ioptions):
    try:
        ok_use = ioptions['cube_orbitals'] or ioptions['comp_p_h_dens'] or ioptions['comp_rho0n']
    except TypeError:
        return ''

    rstr = ''
    if ok_use:
        rstr += addlinec()
        rstr += addlinel("orbkit for orbital/density plotting (http://orbkit.github.io):", 3)
        rstr += addlinel("G. Hermann, V. Pohl, J. C. Tremblay, B. Paulus, H.-C. Hege, A. Schild")
        rstr += addlinel("J. Comput. Chem. (2016), 37, 1511.")

    return rstr

def add_VIST(cfile):
    if cfile is None: return ''

    rstr = ''
    if 'vist' in cfile or 'nics' in cfile:
        rstr += addlinec()
        rstr += addlinel("Visualization of chemical shielding tensors (VIST):", 3)
        rstr += addlinel("F. Plasser, F. Gloecklhofer")
        rstr += addlinel("Eur. J. Org. Chem. (2021), DOI: 10.1002/ejoc.202100352.")

    return rstr
