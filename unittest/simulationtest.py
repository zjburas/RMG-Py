#!/usr/bin/python
# -*- coding: utf-8 -*-
 
import unittest

import sys
sys.path.append('../source')

import pylab
import numpy
import math
		
from rmg.species import *
from rmg.reaction import *
from rmg.model import *

################################################################################

def initializeSimulation1():
	
	speciesA = Species(1, 'A')
	speciesA.thermoData = ThermoGAData(500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
	
	speciesB = Species(2, 'B')
	speciesB.thermoData = ThermoGAData(494237.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
	
	reactionAB = Reaction([speciesA], [speciesB])
	reactionAB.kinetics = [ArrheniusEPKinetics(1.0, 0.0, 0.0, 0.0)]
	
	model = CoreEdgeReactionModel()
	model.addSpeciesToCore(speciesA)
	model.addSpeciesToCore(speciesB)
	model.addReactionToCore(reactionAB)
	model.termination.append(TerminationTime(10.0))
	model.absoluteTolerance = 1e-24
	model.relativeTolerance = 1e-12
	
	system = BatchReactor()
	system.equationOfState = IdealGas()
	system.temperatureModel = TemperatureModel()
	system.temperatureModel.setIsothermal(1000.0)
	system.pressureModel = PressureModel()
	system.pressureModel.setIsobaric(1.0)
	
	system.initialConcentration[speciesA] = 1.0
	
	return speciesA, speciesB, reactionAB, model, system

################################################################################

def simulate(model, system):

	t, y, valid, species = system.simulate(model)
		
	# Reshape y into a matrix rather than a list of lists
	y0 = numpy.zeros((len(t), len(y[0])), float)
	for i, u in enumerate(y):
		for j, v in enumerate(u):
			y0[i,j] = v
			
	return t, y0
	
def postprocess(t, y):

	# Make concentration plot and show
	pylab.plot(t[1:], y[1:,3:5])
	pylab.xlabel('Time (s)')
	pylab.ylabel('Concentration (mol/m^3)')
	pylab.legend(['A', 'B'])
	pylab.show()
	
	print t, y, valid, species

################################################################################

class SimulationCheck(unittest.TestCase):                          

	def testSimulation1(self):
		"""
		Simulation one is a simple isomerization reaction A --> B, with the
		thermodynamics designed for an equilibrium of all B. This occurs in an
		isothermal, isobaric, homogeneous batch reactor.
		"""
	
		speciesA, speciesB, reactionAB, model, system = initializeSimulation1()
		
		# Case 1: Irreversibility
		speciesA.thermoData = ThermoGAData(500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		speciesB.thermoData = ThermoGAData(-500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		
		t, y = simulate(model, system)
		
		# Check equilibrium
		self.assertTrue(y[-1,3] < 0.0001 * y[0,3])
		self.assertTrue(y[-1,4] > 0.9999 * y[0,3])
		
		# Check kinetics
		for i in range(len(t)):
			if abs(t[i] - 1.0) < 0.0001:
				self.assertAlmostEqual(y[i,3] / y[0,3], math.exp(-1.0), 4)
		
	def testSimulation2(self):
		"""
		Simulation two is a simple isomerization reaction A --> B, with the
		thermodynamics designed for an equimolar equilibrium. This occurs in an
		isothermal, isobaric, homogeneous batch reactor.
		"""
	
		speciesA, speciesB, reactionAB, model, system = initializeSimulation1()
		
		# Case 2: Equimolar equilibrium
		speciesA.thermoData = ThermoGAData(500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		speciesB.thermoData = ThermoGAData(500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		
		t, y = simulate(model, system)
			
		# Check equilibrium
		self.assertAlmostEqual(y[-1,3], y[-1,4], 4)
	
	def testSimulation3(self):
		"""
		Simulation three is a simple isomerization reaction A --> B, with the
		thermodynamics designed for an equilibrium ratio of 1:2 A:B. This occurs
		in an isothermal, isobaric, homogeneous batch reactor.
		"""
		speciesA, speciesB, reactionAB, model, system = initializeSimulation1()
	
		# Case 3: A:B = 1:2 at equilibrium
		speciesA.thermoData = ThermoGAData(500000.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		speciesB.thermoData = ThermoGAData(494237.0, 0.0, [0, 0, 0, 0, 0, 0, 0])
		
		t, y = simulate(model, system)
		
		# Check equilibrium
		self.assertAlmostEqual(2.0 * y[-1,3], y[-1,4], 4)
		

		
################################################################################

if __name__ == '__main__':
	unittest.main()