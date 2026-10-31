from .actions import Action
from colt.lazyimport import LazyImportCreator, LazyImporter


with LazyImportCreator() as importer:
    theo_header = importer.lazy_import_as('..theo_header', 'theo_header')
    lib_NICS = importer.lazy_import_as('..lib_NICS', 'lib_NICS')


class NICSBq(Action):

    name = 'nics_bq'

    _colt_description = 'Place ghost atoms (Bq) on a ring axis for NICS calculations'

    _user_input = """
    # Molecular structure in xyz format (Angstrom)
    xyzfile = :: existing_file
    # Ring atoms (starting at 1), in the order they appear around the ring
    ring = :: ilist
    # Heights above the ring centroid (Angstrom), e.g. -H '0 1' for NICS(0) and NICS(1)
    heights = 0 1 :: flist, alias=H
    # Also place points below the ring (negative heights)
    both_sides = false :: bool, alias=b
    # Scan along the ring axis instead: 'hmin hmax step', e.g. -s '0 3 0.5'
    scan = :: flist, optional, alias=s
    # Write the Bq lines to this file (in addition to the screen)
    ofile = :: str, optional, alias=o
    """

    _lazy_imports = LazyImporter({
            '..theo_header': 'theo_header',
            '..lib_NICS': 'lib_NICS',
    })

    def run(xyzfile, ring, heights, both_sides, scan, ofile):
        theo_header.print_header(__class__._colt_description)

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
        print("\nAppend these lines to the molecule specification of the Gaussian input.")
        print("The NICS values are then obtained from the Bq shielding tensors,")
        print("e.g. with: theodore plot_vist -r '%s' <logfile>\n"%' '.join(str(iat) for iat in ring))

        lines = lib_NICS.gaussian_bq_lines(points)
        for line in lines:
            print(line)

        if ofile is not None:
            with open(ofile, 'w') as f:
                f.write('\n'.join(lines) + '\n')
            print("\nBq lines written to %s"%ofile)


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
