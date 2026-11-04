from .actions import Action
from .theotools import timeit
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    lib_NICS = importer.lazy_import_as('..lib_NICS', 'lib_NICS')
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')
    cclib_interface = importer.lazy_import_as('..cclib_interface', 'cclib_interface')
    input_options = importer.lazy_import_as('..input_options', 'input_options')


class PlotVist(Action):

    name = 'plot_vist'

    _user_input = """

    # VIST for only these dummy atoms, e.g. -v '0 3 5'
    vist = :: ilist, optional, alias=v
    # Name of output file (for VMD)
    ofile    = VIST.vmd :: str, alias=o
    # Scale factor for VIST dumb-bells
    scale    = 1.0 :: float, alias=s
    # Create coordinate files (using cclib)
    coor  = false :: bool, alias=c
    # Render and plot all tensors separately
    plot_all = false :: bool, alias=p
    # Add labels for eigenvalues above this value (in ppm), e.g. -l 10
    lab_min = 1000. :: float, alias=l
    # Log files to be parsed
    logfiles = :: list(str)
    # Quantum chemistry code (gaussian, qchem, turbomole, dalton)
    qccode = gaussian :: str, alias=q
    # Ring atoms for NICS(0)/NICS(1)/NICS_zz analysis, e.g. -r '4 5 7 9 8 6' (Gaussian only)
    ring = :: ilist, optional, alias=r
    """

    _colt_description = "Read NICS values and prepare VIST plot"

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..lib_NICS': 'lib_NICS',
            '..error_handler': 'error_handler',
            '..cclib_interface': 'cclib_interface',
            '..input_options': 'input_options',
    })

    @timeit
    def run(vist, ofile, scale, coor, plot_all, lab_min, logfiles, qccode, ring):
        if vist is not None and len(vist) == 0:
            vist = None
        theo_header.print_header('Read NICS values and prepare VIST plot', cfile='plot_VIST.py')

        with open(ofile, 'w') as fh:
            pass

        ioptions = input_options.dens_ana_options(ifile=None, check_init=False)
        ioptions['rtype'] = 'cclib'

        if qccode in ['gaussian', 'g09']: # g09: previous name, still accepted
            nv = lib_NICS.NICS_parser_gaussian()
        elif qccode == 'qchem':
            nv = lib_NICS.NICS_parser_QC()
        elif qccode == 'turbomole':
            nv = lib_NICS.NICS_parser_TM()
        elif qccode == 'dalton':
            nv = lib_NICS.NICS_parser_dalton()
        else:
            raise error_handler.ElseError(qccode, 'qccode')

        for ilog, logfile in enumerate(logfiles):
            nv.read(logfile)
            nv.print_data()
            if ring is not None and len(ring) > 0:
                rfile = 'NICS_ring.txt' if len(logfiles) == 1 else 'NICS_ring%i.txt'%ilog
                nv.ring_analysis(ring, outfile=rfile)
            if coor: # create coor file to be read by VMD
                ioptions['rfile'] = logfile
                ccparser = cclib_interface.file_parser_cclib(ioptions)
                struc = cclib_interface.structure_cclib()
                struc.read_cclib(ccparser.data)
                coorf = "coor%i.xyz"%ilog
                struc.make_coord_file(file_path=coorf,file_type='Bqxyz')
                open(ofile, 'a').write("mol new %s\n"%coorf)
            nv.vmd_tensors(ofile, vist, scale, plot_all, lab_min)

        # Instructions for VMD
        if coor:
            print("""
        Input file for VMD and coordinate file(s) created. Now run:
           vmd -e %s
            """%ofile)
        else:
            print("""
        Input file for VMD created. Now do the following:
        1. Open VMD and load coordinate file
        2.   File - Load Visualization State - %s
            """%ofile)
