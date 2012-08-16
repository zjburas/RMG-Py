#!/usr/bin/python
# -*- coding: utf-8 -*-
################################################################################
#
#   RMG - Reaction Mechanism Generator
#
#   Copyright (c) 2002-2012 Prof. Richard H. West (r.west@neu.edu),
#                           Prof. William H. Green (whgreen@mit.edu)
#                           and the RMG Team (rmg_dev@mit.edu)
#
#   Permission is hereby granted, free of charge, to any person obtaining a
#   copy of this software and associated documentation files (the 'Software'),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be included in
#   all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.
#
################################################################################
import os

import logging

import rmgpy
import rmgpy.qm
import rmgpy.qm.mopac

class QMSettings():
    """
    A minimal class to store settings related to quantum mechanics calculations.
    
    =================== ======================= ====================================
    Attribute           Type                    Description
    =================== ======================= ====================================
    `fileStore`         ``str``                 The path to the QMfiles directory
    `scratchDirectory`  ``str``                 The path to the scratch directory
    =================== ======================= ====================================
    
    """
    def __init__(self):
        self.software = None
        self.fileStore = None
        self.scratchDirectory = None
        self.onlyCyclics = None
        
        RMGpy_path = os.getenv('RMGpy') or os.path.normpath(os.path.join(rmgpy.getPath(),'..'))
        self.RMG_bin_path = os.path.join(RMGpy_path, 'bin')
    
    def checkAllSet(self):
        """
        Check that all the required settings are set.
        """
        assert self.fileStore
        assert self.scratchDirectory
        assert self.software
        assert self.onlyCyclics is not None # but it can be False

class QMCalculator():
    """
    A Quantum Mechanics calculator obect, to store settings. 
    
    The attributes are:

    =================== ======================= ====================================
    Attribute           Type                    Description
    =================== ======================= ====================================
    `settings`          ``QMSettings``          Settings for QM calculations
    =================== ======================= ====================================

    """
    
    def __init__(self,
                 fileStore = None,
                 scratchDirectory = None,
                 ):
        self.settings = QMSettings()
        self.settings.fileStore = fileStore
        self.settings.scratchDirectory = scratchDirectory
        
    def setDefaultOutputDirectory(self, outputDirectory):
        """
        IF the fileStore or scratchDirectory are not already set, put them in here.
        """
        if not self.settings.fileStore:
            self.settings.fileStore = os.path.join(outputDirectory, 'QMfiles')
            logging.info("Setting the quantum mechanics fileStore to {0}".format(self.settings.fileStore))
        if not self.settings.scratchDirectory:
            self.settings.scratchDirectory = os.path.join(outputDirectory, 'QMscratch')
            logging.info("Setting the quantum mechanics scratchDirectory to {0}".format(self.settings.scratchDirectory))
    
    def initialize(self):
        """
        Do any startup tasks.
        """
        self.checkReady()

    def checkReady(self):
        """
        Check that it's ready to run calculations.
        """
        self.settings.checkAllSet()
        self.checkPaths()

    def checkPaths(self):
        """
        Check the paths in the settings are OK. Make folders as necessary.
        """
        self.settings.fileStore = os.path.expandvars(self.settings.fileStore) # to allow things like $HOME or $RMGpy
        self.settings.scratchDirectory = os.path.expandvars(self.settings.scratchDirectory)
        for path in [self.settings.fileStore, self.settings.scratchDirectory]:
            if not os.path.exists(path):
                logging.info("Creating directory %s for QM files."%os.path.abspath(path))
                os.makedirs(path)

        if not os.path.exists(self.settings.RMG_bin_path):
            raise Exception("RMG-Py 'bin' directory {0} does not exist.".format(self.settings.RMG_bin_path))
        if not os.path.isdir(self.settings.RMG_bin_path):
            raise Exception("RMG-Py 'bin' directory {0} is not a directory.".format(self.settings.RMG_bin_path))
            
        
    def getThermoData(self, species):
        """
        Generate thermo data for the given :class:`Species` via a quantum mechanics calculation.
        """
                
        molecule = species.molecule[0]
        
        if self.settings.onlyCyclics and not molecule.isCyclic():
            logging.info("Skipping QM calculation for noncyclic species {0}".format(species.label))
            return None
        
        if self.settings.software == 'mopac':
            qm_molecule_calculator = rmgpy.qm.mopac.MopacMolPM3(molecule, self.settings)
            thermo0 = qm_molecule_calculator.generateThermoData()
        else:
            raise Exception("Unknown QM software '{0}'".format(self.settings.software))
        return thermo0
    
        