Analysis of chemical shielding
------------------------------

TheoDORE provides functionalities to analyse chemical shielding tensors and, in particular their use to analyze aromaticity.


Visualisation of chemical shielding tensors (VIST)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The VIST method [`Eur. J. Org. Chem., (2021) <https://doi.org/10.1002/ejoc.202100352>`_] allows to visualise chemical shielding tensors in the context of the molecular structure showing local variations of the shielding along with its anisotropy. The method proceeds in the following order

* Computation of shielding tensors
* Interfaces exists for Gaussian, Q-Chem,  Turbomole, and Dalton.
* Creation of VIST input via ``theodore plot_vist``
* Rendering of VIST plot in `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_

The tools have only been recently developed and still might require a bit of manual tuning of script files by the user. Any feedback on the user-friendliness of the tools is appreciated.

Command line options
____________________

Run ``theodore plot_vist -h`` to see all command line options

.. code-block:: text

    plot_VIST.py [options] <logfile1> <logfile2>
    Command line options:
    -h/--help            | show this help message and exit                   
    -v/--vist            | VIST for only these dummy atoms, e.g. -v '0 3 5'  
    -o/--ofile           | Name of output file (for VMD)                     
    -s/--scale           | Scale factor for VIST dumb-bells                  
    -c/--coor            | Create coordinate files (using cclib)
    -p/--plot_all        | Render and plot all tensors separately
    -l/--lab_min         | Add labels for eigenvalues above this value (in ppm), e.g. -l 10                                    
    -q/--qccode          | Quantum chemistry code (g09, g16 or gaussian, qchem, turbomole, dalton)
    -r/--ring            | Ring atoms for NICS(0)/NICS(1)/NICS_zz analysis, e.g. -r '4 5 7 9 8 6' (Gaussian only)

Simple example
______________

Create a Gaussian input file for an NMR calculation with the positions of the shielding tensors indicated as dummy (Bq) atoms and run Gaussian.
::

    # PBE1PBE/Def2SVP scf=tight nmr

    Shielding tensors

    0 1
    C          0.00000        1.20802        0.69749
    C          0.00000        1.20802       -0.69749
    C          0.00000       -0.00000       -1.39501
    C         -0.00000       -1.20802       -0.69749
    C         -0.00000       -1.20802        0.69749
    C         -0.00000       -0.00000        1.39501
    H          0.00000        2.15464       -1.24395
    H         -0.00000        2.15464        1.24395
    H         -0.00000       -2.15464       -1.24395
    H          0.00000       -0.00000       -2.48803
    H         -0.00000       -0.00000        2.48803
    H         -0.00000       -2.15464        1.24395
    Bq         0.0            0.0            0.0
    Bq         1.0            0.0            0.0
    Bq         2.0            0.0            0.0


After the job finishes simply run (assuming you named the output file ``gaussian.log``)

::

    plot_VIST.py gaussian.log

Open a molecular coordinate file in VMD. Then from inside VMD click: `File` - `Load Visualization State` - `VIST.vmd`

In the first instance, this will plot the shielding tensors of all Bq positions on top of each other.

.. figure:: figures/VIST-all.png

Use, the `-p` option to plot them all individually or use the `-v` option to select individual tensors.

Processing several output files
_______________________________

In many cases, it is desirable to process several output files at the same time to produce a consistent representation of several output files, for example to represent the shielding in different electronic states. The following figures represent the shielding in the dication, neutral state, and dianion of a macrocycle.


.. image:: figures/2P.png
.. image:: figures/neut.png
.. image:: figures/2M.png




This figure was created with the command

::

    theodore plot_vist -c -v '0 4' */gaussian.log

The ``-c`` option specifies that the molecular coordinates are directly parsed from the log-file. ``-v 0 4`` indicates that only the Bq atoms with index 0 and 4 are included in the representation. ``*/gaussian.log`` indicates that the gaussian.log files in all subdirectories are processed.

The files are loaded into VMD via

::

    vmd -e VIST.vmd

Initially all geometries are superposed on top of each other. To hide some of the molecules use the ``D`` option (part of ``T A D F Molecule``) in the VMD Main menu.

Gaussian: orientation of the ghost atoms
________________________________________

Gaussian reports the shielding tensors in the *standard orientation* unless ``NoSymm`` is given.
``plot_vist`` therefore takes the Bq positions from the last orientation table printed before the
shielding tensors (the standard orientation if present, otherwise the input orientation), so that
coordinates and tensors are always in the same frame.
In practice Gaussian usually turns symmetry off when ghost atoms are present
("Cannot cope with ghost atoms") and prints only the input orientation; both frames then coincide.
The Bq atoms are matched to their tensors through the centre number. This also covers Bq atoms
given as Z-matrix lines or with freeze codes (``Bq 0 x y z``).
If a log file contains several shielding sections (e.g. ``Opt NMR`` or ``--Link1--`` jobs),
only the last one is used. The same parser is used for Gaussian 09 and Gaussian 16.

NICS(0), NICS(1) and NICS_zz
____________________________

Placing the ghost atoms on the ring axis is done with ``theodore nics_bq``.
It reads an xyz file and the ring atoms (in their order around the ring) and prints the Bq lines
to be appended to the Gaussian molecule specification::

    theodore nics_bq molecule.xyz '1 2 3 4 5 6' -H '0 1' -b
    theodore nics_bq molecule.xyz '1 2 3 4 5 6' -s '0 3 0.5'    # scan along the ring axis

The ring centroid is the mean of the ring atom positions and the ring normal is the normal of the
least-squares plane through them, oriented by the right-hand rule with respect to the atom order.

After the NMR job, the ``-r`` option of ``plot_vist`` evaluates all ghost atoms that lie on the
ring axis (within 0.05 A)::

    theodore plot_vist -r '1 2 3 4 5 6' gaussian.log

The table (also written to ``NICS_ring.txt``) lists the height above the ring centroid, the
isotropic value ``NICS_iso`` and the out-of-plane component ``NICS_zz = -n.sigma.n`` with ``n``
the ring normal, and labels the NICS(0) and NICS(+-1) points. NICS values follow the usual
convention NICS = -sigma. The ring atom numbers refer to the order of the real atoms in the log file.
