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

fiveMemberedRingInstance = objects.fiveMemberedRing()

testXYZ = XYZobjects[0]
testXYZ.constructRingIntermediate(fiveMemberedRingInstance)

testXYZ.printToFile("ringTest.xyz")
#testXYZ.viewInVMD()

