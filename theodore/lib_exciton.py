from __future__ import print_function, division

from . import lib_struc, error_handler, units
import numpy

class exciton_analysis:
    """
    Perform analysis of an effective exciton wavefunction.
    Approximate atom centered solutions.
    """
    # TODO: add correlation coefficient and covariance
    #   use the same moment construction as in the JCC paper, only discretized
    
    def __init__(self):
        self.distmat = None
        self.Ohnomat = None
    
    def get_distance_matrix(self, struc, U = 8.0, eps = 2):
        """
        Get distance matrix and compute ERIs according to Ohno formula.
        U ... Hubbard parameter (eV)
        eps ... screening
        """
        self.distmat = struc.ret_distance_matrix()

        # Ohno formula from W. Barford, PRB, 106, 035201 (2022)
        red_distmat = self.distmat / units.energy['eV'] / units.length['A']
        self.Ohnomat = U / (1 + (U * eps * red_distmat)**2)**.5

    def ret_RMSeh(self, Om, OmAt):
        """
        Return the root mean square electron-hole distance (Ang).
        """
        if not type(self.distmat) is numpy.ndarray:
            raise error_handler.MsgError("Compute the distance matrix first!")
        
        try:
            MS_dist = numpy.dot(OmAt.flatten(), self.distmat.flatten()**2.) / Om
        except:
            print("\n Error when computing MS_dist!")
            print(" Please check the coordinate file.")
            print(" OmAt: %i x %i"%(len(OmAt), len(OmAt[0])))
            print(" distmat: %i x %i"%(len(self.distmat), len(self.distmat[0])))
            raise
        
        RMS_dist = numpy.sqrt(MS_dist)
        
        return RMS_dist

    def ret_MAeh(self, Om, OmAt):
        """
        Return the mean absolute electron-hole distance (Ang).
        """
        if not type(self.distmat) is numpy.ndarray:
            raise error_handler.MsgError("Compute the distance matrix first!")
        
        MA_dist = numpy.dot(OmAt.flatten(), self.distmat.flatten()) / Om
        
        return MA_dist

    def ret_K2(self, OmAt):
        """
        Return an approximate K2 exciton binding energy (eV) using the Ohno matrix.
        """
        if not type(self.Ohnomat) is numpy.ndarray:
            raise error_handler.MsgError("Compute the Ohno matrix first!")

        return -numpy.sum(OmAt * self.Ohnomat)

    def ret_rTD(self, QT2, tpop):
        """
        Compute effective transition density size.
        """
        if QT2 < 0.001:
            return None
        pot = 0
        for A, popA in enumerate(tpop):
            for B, popB in enumerate(tpop):
                if A != B:
                    pot += popA * popB / self.distmat[A, B]

        return -QT2 / pot

    def ret_Vinter(self, QT2, tpop):
        """
        Compute interaction term V_inter using Ohno matrix elements.
        """
        if QT2 < 0.001:
            return None
        pot = 0
        for A, popA in enumerate(tpop):
            for B, popB in enumerate(tpop):
                if A != B:
                    pot += popA * popB * self.Ohnomat[A, B]

        return pot

    def ret_Vdiag(self, QT2, tpop):
        """
        Compute diagonal contribution using Ohno matrix.
        """
        if QT2 < 0.001:
            return None
        pot = 0
        for A, popA in enumerate(tpop):
            pot += (popA**2) * self.Ohnomat[A, A]

        return pot

    def ret_J2(self, QT2, tpop):
        """
        Compute self-repulsion of the transition density J2 = V_diag + V_inter.
        - For singlets: real Coulomb self-repulsion.
        - For triplets: formal J2 only (not physically entitled).
        """
        Vdiag = self.ret_Vdiag(QT2, tpop)
        Vinter = self.ret_Vinter(QT2, tpop)

        if Vdiag is None or Vinter is None:
            return None

        J2 = Vdiag + Vinter

        return J2


