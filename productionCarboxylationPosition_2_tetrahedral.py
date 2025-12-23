import functions
import objects

import json
import copy
from tqdm import tqdm
import warnings

#################################################################
#Parameters
#################################################################

fileID = '/home/alpal/projects/methanCapture/carboxylationProblem/db/simah.db'

XYZfile = '/home/alpal/projects/methanCapture/carboxylationProblem/filtered_structures_6256.xyz'

fmaxList = [0.5, 0.25, 0.1, 0.05, 0.025, 0.01]

debug = True

switchR1andR2 = True

outputDir = 'outputs_αR2_tetrahedral'

planar = False

#################################################################
#Main
#################################################################

_, structureCodes = functions.getXYZandBuildingBlock_structureList(fileID, XYZfile, sanityCheck = False)

XYZobjects = functions.load_xyz_structures_from_file(XYZfile)

def xyzObjectCarboxylationRoutine(xyzObject, switchR1andR2 = False, planar = True, outputDir = 'outputs_αR1_planar'):

	fiveMemberedRingInstance = objects.fiveMemberedRing()

	xyzObject.constructRingIntermediate(fiveMemberedRingInstance, switchR1andR2 = False)
	
	if planar:

		xyzObject.forcePlanar()
	else:
		xyzObject.forceTetrahedral()

	xyzObject.reduceRclashes(R_index = 1, rotationAngle = 20)
	xyzObject.reduceRclashes(R_index = 2, rotationAngle = 20)

	xyzObject.pivotCorrectionForValenceSanity(rotationAngle = 20)

	xyzObject.reduceRclashes(R_index = 1, rotationAngle = 20)
	xyzObject.reduceRclashes(R_index = 2, rotationAngle = 20)

	xyzObject.printToFile(f"{outputDir}/ligands_{str(i).zfill(4)}_{structureCodes[i]}.xyz")
	xyzObject.writeSanityRecord(f"{outputDir}/ligands_{str(i).zfill(4)}_{structureCodes[i]}_sanityRecord.txt")

for i, xyzObject in tqdm(enumerate(XYZobjects)):

	xyzObjectCopy = copy.deepcopy(xyzObject)
	converged = False
	fIndex = 0

	try:
		xyzObjectCarboxylationRoutine(xyzObject, switchR1andR2, planar, outputDir)
	except:

		warnings.warn("InitialCarboxlyation failed, trying again with a DFT approximation to fix failed bonds")

		while not converged:

			if fIndex > len(fmaxList):
				print(f"Optimization failed, take note of structure {i}")
				#NOTE: it is "Converged"
				converged = True
				continue

			#NOTE: xyzObjectCopy is not going to be redefined
				#NOTE: ; therefore, it is the original structure from the xyz

			xyzObject = copy.deepcopy(xyzObjectCopy)

			try:

				xyzObject.geometryOptimization(fmax = fmaxList[fIndex])
				xyzObjectCarboxylationRoutine(xyzObject, switchR1andR2, planar, outputDir)
				converged = True

			except:
				print(("carboxylation with optimization failed, trying again with\n"
					"stricter convergence criterion"
				))
				fIndex+=1
				#NOTE: does nothing, but makes my code more readable
				converged = False



