from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    lib_NICS = importer.lazy_import_as('..lib_NICS', 'lib_NICS')
    error_handler = importer.lazy_import_as('..error_handler', 'error_handler')


class NICSPrep(Action):

    name = 'prep_nics'

    _colt_description = 'Place ghost atoms on a ring for NICS/VIST calculations'

    _user_input = """
    # Molecular structure in xyz format (Angstrom)
    xyzfile = :: existing_file
    # Ring atoms (starting at 1), in the order they     appear around the ring, e.g. '1 4 5 6'
    ring = :: ilist
    # Heights above the ring centroid (Angstrom), e.g.  -H '0 1' for NICS(0) and NICS(1)
    heights = 0 1 :: flist, alias=H
    # Also place points below the ring (negative heights)
    both_sides = false :: bool, alias=b
    # Scan along an axis perpendicular to the ring:     'hmin hmax step', e.g. -s '0 3 0.5'
    scan = :: flist, optional, alias=s
    # Write ghost atom lines to this file (in addition  to the screen)
    ofile = :: str, optional, alias=o
    # Quantum chemistry code (currently only gaussian)
    qccode = gaussian :: str, alias=q
    """

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..lib_NICS': 'lib_NICS',
            '..error_handler': 'error_handler',
    })

    def run(xyzfile, ring, heights, both_sides, scan, ofile, qccode):
        theo_header.print_header(__class__._colt_description)

        # General part: ghost atom positions on the ring axis
        coords = read_xyz(xyzfile)
        for iat in ring:
            if not 1 <= iat <= len(coords):
                raise ValueError("Ring atom %i not in %s (%i atoms)"%(iat, xyzfile, len(coords)))
        ring_coords = [coords[iat-1] for iat in ring]

        if scan is not None and len(scan) > 0:
            if len(scan) != 3:
                raise ValueError("scan needs three values: hmin hmax step")
            points, hlist = lib_NICS.ring_scan_points(ring_coords, *scan)
        else:
            points, hlist = lib_NICS.ring_probe_points(ring_coords, heights, both_sides)

        centroid, normal = lib_NICS.ring_frame(ring_coords)
        print("Ring centroid: (% .5f, % .5f, % .5f)"%tuple(centroid))
        print("Ring normal:   (% .5f, % .5f, % .5f)"%tuple(normal))
        print("Heights (A):", ' '.join('%.3f'%h for h in hlist))

        # Program specific part: input lines for the ghost atoms
        if qccode in ['gaussian', 'g09']:
            lines = lib_NICS.gaussian_bq_lines(points)
            print("\nAppend these lines to the molecule specification of the Gaussian input.")
        else:
            raise error_handler.ElseError(qccode, 'qccode')

        print("The NICS values are then obtained from the shielding tensors of the ghost atoms,")
        print("e.g. with: theodore plot_vist -q %s -r '%s' <logfile>\n"%(qccode, ' '.join(str(iat) for iat in ring)))

        for line in lines:
            print(line)

        if ofile is not None:
            with open(ofile, 'w') as f:
                f.write('\n'.join(lines) + '\n')
            print("\nGhost atom lines written to %s"%ofile)


def read_xyz(xyzfile):
    """
    Read the coordinates (Angstrom) of the first structure in an xyz file.
    """
    with open(xyzfile, 'r') as f:
        natom = int(f.readline().split()[0])
        f.readline()
        coords = []
        for iat in range(natom):
            words = f.readline().split()
            coords.append([float(w) for w in words[1:4]])
    return coords
