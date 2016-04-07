# Data sources		
database(		
    thermoLibraries = ['KlippensteinH2O2','primaryThermoLibrary','DFT_QCI_thermo','CBS_QB3_1dHR','Chernov'],		
    reactionLibraries = [('Chernov',False)],		
    seedMechanisms = ['KlippensteinH2O2','ERC-FoundationFuelv0.9'],		
    kineticsDepositories = ['training'],		
    kineticsFamilies = 'default',		
    kineticsEstimator = 'rate rules',		
)		
		
generatedSpeciesConstraints(		
    allowed=['seed mechanisms','input species','reaction libraries'],		
    maximumOxygenAtoms=4,	
    maximumCarbonAtoms=6,		
    maximumRadicalElectrons=2,	
    allowSingletO2 = True,	
)		
		
		
# List of species		
species(		
    label='CH4',		
    reactive=True,		
    structure=SMILES("C"),		
)		
		
species(		
    label='O2',		
    reactive=True,		
    structure=SMILES('[O][O]'),		
)		
		
species(		
    label='N2',		
    reactive=False,		
    structure=SMILES("N#N"),		
)		
		
species(		
    label='AC3H4',		
    reactive=True,		
        structure = adjacencyList(		
  """		
  1 C u0 p0 c0 {2,D} {4,S} {5,S}		
  2 C u0 p0 c0 {1,D} {3,D}		
  3 C u0 p0 c0 {2,D} {6,S} {7,S}		
  4 H u0 p0 c0 {1,S}		
  5 H u0 p0 c0 {1,S}		
  6 H u0 p0 c0 {3,S}		
  7 H u0 p0 c0 {3,S}		
  """),		
)		
		
		
species(		
    label='C3H8',		
    reactive=True,		
        structure = adjacencyList(		
  """		
  1  C u0 p0 c0 {2,S} {4,S} {5,S} {6,S}		
  2  C u0 p0 c0 {1,S} {3,S} {10,S} {11,S}		
  3  C u0 p0 c0 {2,S} {7,S} {8,S} {9,S}		
  4  H u0 p0 c0 {1,S}		
  5  H u0 p0 c0 {1,S}		
  6  H u0 p0 c0 {1,S}		
  7  H u0 p0 c0 {3,S}		
  8  H u0 p0 c0 {3,S}		
  9  H u0 p0 c0 {3,S}		
  10 H u0 p0 c0 {2,S}		
  11 H u0 p0 c0 {2,S}		
  """),		
)		
		
species(		
    label='CO2',		
    reactive=True,		
    structure=SMILES('O=C=O'),		
)		
		
species(		
    label='C2H6',		
    reactive=True,		
    structure=SMILES('CC'),		
)		
		
species(		
    label='C4H10',		
    reactive=True,		
    structure=SMILES('CCCC'),		
)		
		
species(		
    label='iC5H12',		
    reactive=True,		
    structure=SMILES('CCC(C)C'),		
)		
		
species(		
    label='C6H14',		
    reactive=True,		
    structure=SMILES('CCCCCC'),		
)		
		
		
# Reaction systems		
		
simpleReactor(		
    temperature=(800,'K'),		
    pressure=(2.5,'atm'),		
    initialMoleFractions={		
   "CH4": 	0.343116861	,
   "O2":	0.63	,
   "AC3H4":	0.000133974	,
   "C3H8":	0.000996137	,
   "CO2": 	0.005203451	,
   "C2H6": 	0.018870629	,
   "C4H10": 	0.000100633	,
   "C6H14":	5.74419E-05	,
   "iC5H12": 	0.0000207	,
   "N2":	0.006009077	,
},		
terminationTime=(0.05,'s'),		
)

simpleReactor(		
    temperature=(1200,'K'),		
    pressure=(2.5,'atm'),		
    initialMoleFractions={		
   "CH4": 	0.343116861	,
   "O2":	0.63	,
   "AC3H4":	0.000133974	,
   "C3H8":	0.000996137	,
   "CO2": 	0.005203451	,
   "C2H6": 	0.018870629	,
   "C4H10": 	0.000100633	,
   "C6H14":	5.74419E-05	,
   "iC5H12": 	0.0000207	,
   "N2":	0.006009077	,
},		
terminationTime=(0.05,'s'),		
)

simpleReactor(		
    temperature=(1600,'K'),		
    pressure=(2.5,'atm'),		
    initialMoleFractions={		
   "CH4": 	0.343116861	,
   "O2":	0.63	,
   "AC3H4":	0.000133974	,
   "C3H8":	0.000996137	,
   "CO2": 	0.005203451	,
   "C2H6": 	0.018870629	,
   "C4H10": 	0.000100633	,
   "C6H14":	5.74419E-05	,
   "iC5H12": 	0.0000207	,
   "N2":	0.006009077	,
},		
terminationTime=(0.05,'s'),		
)

simpleReactor(		
    temperature=(2000,'K'),		
    pressure=(2.5,'atm'),		
    initialMoleFractions={		
   "CH4": 	0.343116861	,
   "O2":	0.63	,
   "AC3H4":	0.000133974	,
   "C3H8":	0.000996137	,
   "CO2": 	0.005203451	,
   "C2H6": 	0.018870629	,
   "C4H10": 	0.000100633	,
   "C6H14":	5.74419E-05	,
   "iC5H12": 	0.0000207	,
   "N2":	0.006009077	,
},		
terminationTime=(0.05,'s'),		
)

simpleReactor(		
    temperature=(2400,'K'),		
    pressure=(2.5,'atm'),		
    initialMoleFractions={		
   "CH4": 	0.343116861	,
   "O2":	0.63	,
   "AC3H4":	0.000133974	,
   "C3H8":	0.000996137	,
   "CO2": 	0.005203451	,
   "C2H6": 	0.018870629	,
   "C4H10": 	0.000100633	,
   "C6H14":	5.74419E-05	,
   "iC5H12": 	0.0000207	,
   "N2":	0.006009077	,
},		
terminationTime=(0.05,'s'),		
)

simpleReactor(		
    temperature=(2700,'K'),		
    pressure=(2.5,'atm'),		
    initialMoleFractions={		
   "CH4": 	0.343116861	,
   "O2":	0.63	,
   "AC3H4":	0.000133974	,
   "C3H8":	0.000996137	,
   "CO2": 	0.005203451	,
   "C2H6": 	0.018870629	,
   "C4H10": 	0.000100633	,
   "C6H14":	5.74419E-05	,
   "iC5H12": 	0.0000207	,
   "N2":	0.006009077	,
},		
terminationTime=(0.05,'s'),		
)


model(		
    toleranceKeepInEdge=0.01,		
    toleranceMoveToCore=0.1,		
    toleranceInterruptSimulation=10^8,		
    maximumEdgeSpecies=100000,
    #filterReactions=True,	
)		
		
pressureDependence(		
    method='modified strong collision',		
    maximumGrainSize=(0.5, 'kcal/mol'),		
    minimumNumberOfGrains=250,		
    temperatures=(300,2750,'K',8),		
    pressures=(0.01,100,'bar',5),		
    interpolation=('Chebyshev', 6, 4),		
    maximumAtoms=16,		
)
		
options(		
    units='si',		
    saveRestartPeriod=None,		
    saveSimulationProfiles=True,		
    generateOutputHTML=False,		
    generatePlots=False,		
)		
