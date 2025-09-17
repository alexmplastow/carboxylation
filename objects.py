import functions

import tempfile
import subprocess
import os
import numpy as np



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



		




		
		



