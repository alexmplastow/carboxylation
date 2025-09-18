import functions

import tempfile
import subprocess
import os
import numpy as np

from ase.io import read
from ase.neighborlist import NeighborList, natural_cutoffs

#TODO: add a method for fishing out bond lengths from the example structure you constructed in gaussview
#TODO: but there is no reason you can't put them in their own data structure
#TODO: add a method for distancing C1 and C2 from one another
#TODO: add a method to the carboxylate method for fixing bond lengths
#TODO: add a method for paritionioning R groups from the xyz file

#Add a routine for adding the pentane ring to the molecule of interest

class atom:
	def __init__(self, atomString):
		self.atomString = atomString
		self.element = self.atomString.split(" ")[0]
		self.r = np.array(self.atomString.split(" ")[1:], dtype = np.float32)


class xyzStructure:
	def __init__(self, xyzString):
		self.xyzString = xyzString
		lines = xyzString.split("\n")
		self.atomCount = int(lines[0])
		self.atoms = [ln.strip() for ln in lines[2:] if ln.strip()]

	def printToFile(self, xyzPath):
		file = open(xyzPath, "a")
		file.write(self.xyzString)
		file.close()

	def viewInVMD(self):
		# Courtesy of chatGPT
        	with tempfile.NamedTemporaryFile(delete=False, suffix=".xyz") as tmp:
            		tmp.write(self.xyzString.encode("utf-8"))
            		tmp_path = tmp.name

        	try:
            		subprocess.run(["vmd", tmp_path])

        	finally:
            		if os.path.exists(tmp_path):
                		os.remove(tmp_path)

	def returnAtomByIndex(self, i):
		atomLines = [line for line in self.xyzString.split("\n") 
				if len(line.split(" ")) == 4]
		return atom(atomLines[i])

	def findX1_andX2nearNickelIndices(self, prefer_two_shortest=True, scale=1.25, metal_fudge=0.20, particleType="c"):
		xzy_text = self.xyzString
		atoms = functions.parse_xyz(xzy_text)

		if not atoms:
			raise ValueError("No atoms parsed from XYZ.")

		symbols = [s for s, _ in atoms]
		positions = [p for _, p in atoms]
		adj = functions.build_connectivity(symbols, positions, scale=scale, metal_fudge=metal_fudge)

		ni_list = [i for i, s in enumerate(symbols) if s.lower() == "ni"]

		if not ni_list:
			raise ValueError("No Ni atom found.")

		results = []

		for i_ni in ni_list:
			# carbon neighbors of Ni
			c_neighbors = [j for j in adj[i_ni] if symbols[j].lower() == particleType]
			if len(c_neighbors) < 2:
				results.append({
					"Ni_index": i_ni,
					"error": f"Found only {len(c_neighbors)} carbon neighbor(s) for Ni; need at least 2."
				})
				continue

			# choose the two closest if more than 2
			c_neighbors_sorted = sorted(c_neighbors, key=lambda j: functions.dist(positions[i_ni], positions[j]))

			# NOTE: c_neighbors_sorted does not look empty
			if prefer_two_shortest:
				C1, C2 = c_neighbors_sorted[:2]
			else:
				C1, C2 = c_neighbors_sorted[0], c_neighbors_sorted[1]

			# Ensure consistent ordering by index (optional)
			if C2 < C1:
				C1, C2 = C2, C1

			return C1, C2

#Write a method to pluck a specific atom index from your XYZ files

	#NOTE: maybe rename this one in the event that you switch metal types
	def findX1_andX2nearNickelRs(self, particleType="c"):
		C1_i, C2_i = self.findX1_andX2nearNickelIndices(particleType=particleType)

		C1 = self.returnAtomByIndex(C1_i)
		C2 = self.returnAtomByIndex(C2_i)

		C1_r = C1.r
		C2_r = C2.r

		return C1_r, C2_r
	def findNickelR(self):
		metalAtom_index = [i for i, atomLine in enumerate(self.atoms) 
					if "Ni" in atomLine][0]
		metalAtom = self.returnAtomByIndex(metalAtom_index)
		Ni_r = metalAtom.r
		return Ni_r

	#TODO: add function for attaching the carboxylate to the metal
	def addCarboxylate(self, carboxylate):
		self.carboxylate = carboxylate

		self.atomCount+=3

		self.atoms = 	self.atoms + \
				[carboxylate.O1line,
				carboxylate.O2line, 
				carboxylate.Cline]

		self.xyzString = str(self.atomCount) + "\n" + "\n" + "\n".join(self.atoms)

	#Note: I am not 100% sure this code will generalize well
	def getRGroupIndices(self, index = 1):

		if index == 1:

			C1_i, C2_i = self.findX1_andX2nearNickelIndices()

			R_C = self.returnAtomByIndex(C1_i)
			R_C_i = C1_i
			R_C_j = C2_i

		elif index == 2:

			C1_i, C2_i = self.findX1_andX2nearNickelIndices()

			R_C = self.returnAtomByIndex(C2_i)
			R_C_i = C2_i
			R_C_j = C1_i
		else:
			print("You had to pick 1 or 2 as the input, was that so hard?")
			raise Exception("Bruh")
		
		#TODO: MAKE SURE TO DELETE THIS FILE WHEN YOU ARE DONE
		self.printToFile("tmp.xyz")


		atoms = read("tmp.xyz")
		cutoffs = natural_cutoffs(atoms)
		nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
		nl.update(atoms)	
		indices, _ = nl.get_neighbors(R_C_i)
		indices = [int(index) for index in list(indices)]
		R_indices = [int(index) for index in indices if index != R_C_j and index != 0]
		R_indices.append(R_C_i)
		checkedR_indices = [R_C_i]

		while R_indices != checkedR_indices:
			R_indices = sorted(R_indices)
			checkedR_indices = sorted(checkedR_indices)
			
			for R_index in R_indices:
				if R_index not in checkedR_indices:
					neighborIndices, _ = nl.get_neighbors(R_index)
					
					checkedR_indices.append(int(R_index))

					for neighborIndex in neighborIndices:
						if neighborIndex not in R_indices:
							R_indices.append(int(neighborIndex))

				else:
					continue
		os.remove("tmp.xyz")


		return R_indices

	def getRGroup(self, index = 1, deleteRgroupFromXYZ = True):
		R_indices = self.getRGroupIndices(index = index)
		R_xyz = [self.atoms[i] for i in R_indices]
		if deleteRgroupFromXYZ == True:
			self.atoms = [self.atoms[i] for i in 
					range(0,len(self.atoms)) 
					if i not in R_indices]
			self.atomCount = str(int(self.atomCount) - len(R_indices))
			#Making a change to the xyz string
			self.xyzString = self.atomCount + "\n" + "\n" + "\n".join(self.atoms)
		return Rgroup(R_indices, R_xyz)

			

class Rgroup:
	def __init__(self, R_indices, R_xyz):
		self.R_indices = R_indices
		#Note: this attribute contains some lines from the xyz file
			#Note: it does not contain all the information needed for a .xyz
		self.R_xyz = R_xyz
		
	

#TODO: clean up the redefining of lines a bit, this is what functions are made for
		
		
class carboxylate:
	def __init__(self, xyzPath = "/home/alpal/projects/methanCapture/carboxylationProblem/data/carboxylate.xyz"):
		self.xyzPath = xyzPath
		parsedXYZ = functions.parse_xyz(open(self.xyzPath, "r").read())
		self.O1 = parsedXYZ[0]
		self.O2 = parsedXYZ[1]
		self.C = parsedXYZ[2]
		
		self.r_O1 = np.array(self.O1[-1])
		self.r_O2 = np.array(self.O2[-1])
		self.r_C = np.array(self.C[-1])
		
		self.O1line = f"{self.O1[0]} {self.r_O1[0]} {self.r_O1[1]} {self.r_O1[2]}\n"
		#Note: I changed this line on 9/17/25
		self.O2line = f"{self.O2[0]} {self.r_O2[0]} {self.r_O2[1]} {self.r_O2[2]}\n"
		self.Cline = f"{self.C[0]} {self.r_C[0]} {self.r_C[1]} {self.r_C[2]}\n"


	#Note: courtesy of chatGPT
	def putInPlaneOf_MetalAlkyne(self, xyzObject):
		src_pts = np.vstack([self.r_O1, self.r_O2, self.r_C])

		C1_r, C2_r = xyzObject.findX1_andX2nearNickelRs()
		Ni_r = xyzObject.findNickelR()

		dst_pts = np.vstack([np.array(C1_r), np.array(C2_r), np.array(Ni_r)])

		R, t = functions.rigid_transform(src_pts, dst_pts)

		# Apply
		self.r_O1 = R @ self.r_O1 + t
		self.r_O2 = R @ self.r_O2 + t
		self.r_C  = R @ self.r_C  + t

		# Update lines (if you want to write back to xyz)
		self.O1line = f"{self.O1[0]} {self.r_O1[0]} {self.r_O1[1]} {self.r_O1[2]}\n"
		self.O2line = f"{self.O2[0]} {self.r_O2[0]} {self.r_O2[1]} {self.r_O2[2]}\n"
		self.Cline  = f"{self.C[0]} {self.r_C[0]} {self.r_C[1]} {self.r_C[2]}\n"

	def printToXYZ(self, xyzPath):
		self.xyzString = "3\n\n" + self.O1line + self.O2line + self.Cline
		file = open(xyzPath, "a")
		file.write(self.xyzString)
		file.close()
	
	#NOTE: while I technically do not need the entire molecule just to rotate
		#NOTE: the carboxylate 
	def rotateAboutC(self, xyzObject, angle):

		C1_r, C2_r = xyzObject.findX1_andX2nearNickelRs()
		Ni_r = xyzObject.findNickelR()

		cross = functions.findCrossProduct(Ni_r, C1_r, C2_r)

		rotatedPoints = functions.rotatePointsAboutCrossProduct(
					[self.r_O1, self.r_O2, self.r_C],
					cross,
					angle,
					pivot_index = -1)

		self.r_O1 = rotatedPoints[0]
		self.r_O2 = rotatedPoints[1]
		self.r_C = rotatedPoints[2]

		self.O1line = f"{self.O1[0]} {self.r_O1[0]} {self.r_O1[1]} {self.r_O1[2]}\n"
		self.O2line = f"{self.O2[0]} {self.r_O2[0]} {self.r_O2[1]} {self.r_O2[2]}\n"
		self.Cline  = f"{self.C[0]} {self.r_C[0]} {self.r_C[1]} {self.r_C[2]}\n"


#TODO: add a method to the xyzStructure for identifying the indices of those atoms which are attached to C1 and C2 (this could be a headache..... I'll probably need a molecular graphics software to guess bonds for me

#TODO: add an R1andR2 object for handling rotation operations and such
	#TODO: it is best if you can add these as attributes to the xyzStructure



		
		



