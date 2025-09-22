import functions
import objects

import json

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


for i, xyzObject in enumerate(XYZobjects):
	try:

		fiveMemberedRingInstance = objects.fiveMemberedRing()

		xyzObject.constructRingIntermediate(fiveMemberedRingInstance, switchR1andR2 = True)
		xyzObject.separateFiveMemberRing(fiveMemberedRingInstance, d_sep = 0.5)
		xyzObject.printToFile(f"outputs_αR2/ligands_{str(i).zfill(4)}.xyz")

	except:
		print("I haven't configured most of these routines to handle cases of ligand complexes without one of the Ni branches, I'm assuming that's what this is")

