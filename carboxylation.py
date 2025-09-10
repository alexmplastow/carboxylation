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



