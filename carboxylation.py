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
testXYZ.addCarboxylate(carboxylateInstance)

#testXYZ.viewInVMD()

#testXYZ.add_carboxylate_to_C1()

# Or, alternatively, to C2
# testXYZ.add_carboxylate_to_C2()

# Visualize the result
#testXYZ.viewInVMD()

# If you want to inspect the modified XYZ text:
#print(testXYZ.xyzString)




