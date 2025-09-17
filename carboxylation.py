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

testXYZ = XYZobjects[0]

carboxylateInstance = objects.carboxylate()
carboxylateInstance.putInPlaneOf_MetalAlkyne(testXYZ)
carboxylateInstance.printToXYZ("carboxylateTest.xyz")

carboxylateInstance.rotateAboutC(testXYZ, angle = 30)

carboxylateInstance.printToXYZ("carboxylateTest_rotation.xyz")




