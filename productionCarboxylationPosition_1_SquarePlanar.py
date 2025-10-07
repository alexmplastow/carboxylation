import functions
import objects

import json
import copy

#################################################################
#Parameters
#################################################################

fileID = '/home/alpal/projects/methanCapture/carboxylationProblem/db/simah.db'

#################################################################
#Main
#################################################################

#NOTE: These are not file paths, they are just strings
xyzList = functions.getXYZstructureList(fileID)
XYZobjects = [objects.xyzStructure(xyzString) for xyzString in xyzList]

#xyzObject = XYZobjects[0]

#fiveMemberedRingInstance = objects.fiveMemberedRing()
	
#xyzObjectOrigninal = copy.deepcopy(xyzObject)

#xyzObject.viewInVMD()
#xyzObject.constructRingIntermediate(fiveMemberedRingInstance, switchR1andR2 = False)


for i, xyzObject in enumerate(XYZobjects):

	fiveMemberedRingInstance = objects.fiveMemberedRing()
	
	xyzObjectOrigninal = copy.deepcopy(xyzObject)

	try:
		xyzObject.constructRingIntermediate(fiveMemberedRingInstance, switchR1andR2 = False)
	
	except:

		print("I haven't configured most of these routines to handle cases of ligand complexes without one of the Ni branches, I'm assuming that's what this is")
		continue
	try:
		xyzObject.forcePlanar()
		xyzObject.reduceRclashes(R_index = 1)
		xyzObject.reduceRclashes(R_index = 2)
		xyzObject.pivotCorrectionForValenceSanity()
		xyzObject.printToFile(f"outputs_αR1_planar/ligands_{str(i).zfill(4)}.xyz")
		xyzObject.writeSanityRecord(f"outputs_αR1_planar/ligands_{str(i).zfill(4)}_sanityRecord.txt", hydrogenSkip = True)

	except:
		print("That's not good")
	
	
