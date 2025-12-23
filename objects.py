#NOTE: this code is also rather delicate with regard to radians and degrees, some are in radians other are in degrees. Right now that must be inferred from the code which is acutely lame

#NOTE: the code here is fairly intuitive, spare the fact that atom.r and atom.atomString update 
	#NOTE: seperately

#TODO: work out a fix for enforcing bond lengths closer to your analog (i.e. maybe use an otpimizer to translate the five member ring intermediate until the bond lengths are close to the analog)

#TODO: add a scheme for generating both the tetrahedral and planar configurations in addition to the C1C2Ni plane alignment (best to work Thang's code into your xyzStructure object's methods)
	#NOTE: I think it is fairly clear that I can solve the problem a little faster by working out
	#NOTE: my own transformations instead of relying on Thang's code

#NOTE: at least one of the molecules optimized isn't sane, 16 looks like it has a double bonded hydrogen,
	#NOTE: even before the molecule had a hydrogen added

#NOTE: ligand 250 might be another good example, it seems like the molecule in question has a couple
	#NOTE: close hydrogens, VMD doesn't perceive they are bonds, but ASE does
#TODO: the tetrahedral fix is slightly imperfect, fix it if you find the time

#TODO: your carboxylate vanishes after R rotation, I have no idea why


#TODO: work out a method which handles reorienations of the R1/R2 groups to avoid steric clashes
	#NOTE: 3 is a good test bed for this
	#NOTE:  4 is another good case
	#NOTE: so is 16
	#NOTE: actually this seems to be the issue in pretty much all of these guys

import functions

import tempfile
import subprocess
import shutil
import os
import numpy as np

import math

from ase.io import read
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.optimize import BFGS
from ase.calculators.emt import EMT

from rdkit.Chem.rdmolfiles import MolFromXYZFile
from rdkit.Chem import rdmolops
from rdkit.Chem import rdDetermineBonds
from rdkit import Chem



MAX_VALENCE = {
    "H": 1,
    "He": 0,
    "Li": 1,
    "Be": 2,
    "B": 3,
    "C": 4,
    "N": 3,
    "O": 2,
    "F": 1,
    "Ne": 0,
    "Na": 1,
    "Mg": 2,
    "Al": 3,
    "Si": 4,
    "P": 3,
    "S": 2,
    "Cl": 1,
    "Ar": 0,
    "K": 1,
    "Ca": 2,
    "Sc": 3,
    "Ti": 4,
    "V": 5,
    "Cr": 6,
    "Mn": 7,
    "Fe": 6,
    "Co": 6,
    "Ni": 6,
    "Cu": 4,
    "Zn": 2,
    "Ga": 3,
    "Ge": 4,
    "As": 3,
    "Se": 2,
    "Br": 1,
    "Kr": 0,
    "Rb": 1,
    "Sr": 2,
    "Y": 3,
    "Zr": 4,
    "Nb": 5,
    "Mo": 6,
    "Tc": 7,
    "Ru": 8,
    "Rh": 6,
    "Pd": 6,
    "Ag": 4,
    "Cd": 2,
    "In": 3,
    "Sn": 4,
    "Sb": 3,
    "Te": 2,
    "I": 1,
    "Xe": 0
}



#Add a routine for adding the pentane ring to the molecule of interest

#TODO: generate the molecules

#NOTE:
class atom:
	def __init__(self, atomString):
		self.atomString = atomString
		self.element = self.atomString.split(" ")[0]
		self.r = np.array([part for part in self.atomString.split(" ")[1:] 
					if part != ''], dtype = np.float32)
	#NOTE: changing the r attribute does not directly change the representative string
	def updateString(self):
		self.atomString = f"{self.element} {self.r[0]} {self.r[1]} {self.r[2]}"

	def translate(self, translationVector):
		self.r += translationVector
		self.atomString = f"{self.element} {self.r[0]} {self.r[1]} {self.r[2]}"


class xyzStructure:
	def __init__(self, xyzString):
		self.xyzString = xyzString
		lines = xyzString.split("\n")
		self.atomCount = int(lines[0])
		self.name_line = lines[1].strip() 
		self.atomLines = [ln.strip() for ln in lines[2:] if ln.strip()]
		self.atoms = [atom(atomLine) for atomLine in self.atomLines]

	#XTB is fitted 
	def geometryOptimization(self, fmax = 0.50):

		from xtb.ase.calculator import XTB
		from ase.optimize import BFGS

		tmpFileName = functions.random_filename(prefix='tmp', suffix=".xyz", length=8)
		self.printToFile(tmpFileName)

		aseAtoms = read(tmpFileName)
		aseAtoms.calc = XTB(method="GFN2-xTB")
		optimizer = BFGS(aseAtoms)
		optimizer.run(fmax= fmax)

		optimizedAtoms = []

		for aseAtom in aseAtoms:
			element = aseAtom.symbol
			x, y, z = aseAtom.position
			atomString = f"{element} {x:.8f} {y:.8f} {z:.8f}"
			atomInstance = atom(atomString)
			optimizedAtoms.append(atomInstance)

		self.regenerateAtomLines(optimizedAtoms)

		'''
		tmpFileName = functions.random_filename(prefix='tmp', suffix=".xyz", length=8)
		self.printToFile(tmpFileName)
		#NOTE: this is not an attribute
		aseAtoms = read(tmpFileName)
		aseAtoms.calc = EMT()
		optimizer = BFGS(aseAtoms)
		optimizer.run(fmax = 0.05)

		optimizedAtoms = []

		for aseAtom in aseAtoms:

			element = aseAtom.symbol
			r = aseAtom.position

			atomString = f"{element} {r[0]} {r[1]} {r[2]}"
			atomInstance = atom(atomString)

			optimizedAtoms.append(atomInstance)

		self.regenerateAtomLines(optimizedAtoms)
		'''
	
	#NOTE: assumes the atoms attribute has been set correctly
		#Fixes the xyz string while its at it
	#NOTE: as such, I've forced the atoms attribute to be a manual input
	def regenerateAtomLines(self, atoms):
		
		for atom in self.atoms:
			atom.updateString()

		self.atomCount = len(atoms)
		self.atomLines = [atom.atomString for atom in atoms]
		self.xyzString = str(self.atomCount) + "\n" + "\n" + "\n".join(self.atomLines)


	def printToFile(self, xyzPath):
		file = open(xyzPath, "a")
		file.write(self.xyzString)
		file.close()


	def _view_temp_xyz(self, viewer):
		with tempfile.NamedTemporaryFile(prefix="view_", suffix=".xyz", delete=False) as tmp:
			tmp.write(self.xyzString.encode("utf-8"))
			tmp_path = tmp.name
		try:
			subprocess.run([viewer, tmp_path])
		finally:
			if os.path.exists(tmp_path):
				os.remove(tmp_path)

	def viewInVMD(self):
		self._view_temp_xyz("vmd")

	def viewInAvogadro(self):
		self._view_temp_xyz("avogadro")


	def returnAtomByIndex(self, i):
		atomLines = []
		for line in self.xyzString.split("\n"):
			parts = line.split(" ")
			parts = [part for part in parts if part != '']
			if len(parts) == 4:
				atomLines.append(line)
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
		metalAtom_index = [i for i, atomLine in enumerate(self.atomLines) 
					if "Ni" in atomLine][0]
		metalAtom = self.returnAtomByIndex(metalAtom_index)
		Ni_r = metalAtom.r
		return Ni_r

	#TODO: add function for attaching the carboxylate to the metal
	def addCarboxylate(self, carboxylate):
		self.carboxylate = carboxylate

		self.atomCount+=3

		self.atomLines = self.atomLines + \
				[carboxylate.O1line,
				carboxylate.O2line, 
				carboxylate.Cline]

		self.xyzString = str(self.atomCount) + "\n" + "\n" + "\n".join(self.atomLines)

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
		tmpFileName = functions.random_filename(prefix = 'tmp', suffix='.xyz')
		self.printToFile(tmpFileName)


		atoms = read(tmpFileName)
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
		os.remove(tmpFileName)


		return R_indices

	def getRGroup(self, index = 1, deleteRgroupFromXYZ = True):
		
		C1_r, C2_r = self.findX1_andX2nearNickelRs()
		if index == 1:
			C_r = C1_r
		elif index == 2:
			C_r = C2_r

		R_indices = self.getRGroupIndices(index = index)
		#Note: I changed this from a list, it is probs fine
		R_xyz = "\n".join([self.atomLines[i] for i in R_indices])
		if deleteRgroupFromXYZ == True:
			self.atomLines = [self.atomLines[i] for i in 
					range(0,len(self.atomLines)) 
					if i not in R_indices]
			self.atomCount = str(int(self.atomCount) - len(R_indices))
			#Making a change to the xyz string
			self.xyzString = self.atomCount + "\n" + "\n" + "\n".join(self.atomLines)
		return Rgroup(R_indices, R_xyz, C_r)


	#TODO: extend this to handle cases where the carboxylate has been added
	def deleteRgroups(self, carboxylated = False):
		if not carboxylated:
			R1_indices = self.getRGroupIndices(index = 1)
			R2_indices = self.getRGroupIndices(index = 2)

		elif carboxylated:
			R1_indices = self.R1_indices
			R2_indices = self.R2_indices

		else:
			raise Exception("🤦🤦🤦🤦🤦🤦🤦🤦🤦🤦🤦🤦🤦🤦")

		R_indices = R1_indices + R2_indices




		self.atomLines = [self.atomLines[i] for i in 							range(0,len(self.atomLines)) 
					if i not in R_indices]

		self.atoms = [atom(self.atomLines[i]) for i in 
				range(0,len(self.atomLines)) 
				if i not in R_indices]


		self.atomCount = str(int(self.atomCount) - len(R_indices))
		self.xyzString = self.atomCount + "\n" + "\n" + "\n".join(self.atomLines)


	#TODO: this function still needs tested
	#NOTE: I may not actually need it
	def changeC1C2distance(self, fiveMemberedRing):
		#Find the midpoint between the two carbons

		C1_index, C2_index = self.findX1_andX2nearNickelIndices()
		C1_r, C2_r = self.findX1_andX2nearNickelRs()
		m = self.findAtomByIndex(C1_index, C2_index)

		#Use that to define the direction of their unit vector
		C1_m_uv = functions.unitVector(C1_r, m)
		C2_m_uv = functions.unitVector(C2_r, m)

		#Grab their respective R groups
		R1 = getRGroup(index = 1, deleteRgroupFromXYZ = False)
		R2 = getRGroup(index = 2, deleteRgroupFromXYZ = False)

		#Delete the respectives indices from the xyzString
		self.deleteRgroups()		

		#Grab the new distance, which I am assuming is from the five membered ring
		newD = fiveMemberedRing.d_C1C2

		#Translate both R groups by the (newDistance - oldDistance)/2
		oldD = self.findAtom_AtomDistance(C1_index, C2_index)

		dD = (newD - oldD)/2

		R1.translate(dB * C1_m_uv)
		R2.translate(dB * C2_m_uv)

		#I need to add the R groups back to the xyzString
		self.xyzString += R1.R_xyz
		self.xyzString += R2.R_xyz

		#Redefining the atomLines attribute
		lines = xyzString.split("\n")
		self.atomLines = [ln.strip() for ln in lines[2:] if ln.strip()]
		
		
	def findAtom_AtomDistance(self, atom1Index, atom2Index):
		
		atom1 = self.returnAtomByIndex(atom1Index)
		atom2 = self.returnAtomByIndex(atom2Index)

		d = np.linalg.norm(atom1.r - atom2.r)
		return d

	def findAtom_AtomMidpoint(self, atom1Index, atom2Index):

		atom1 = self.returnAtomByIndex(atom1Index)
		atom2 = self.returnAtomByIndex(atom2Index)

		m = (atom1.r + atom2.r)/2
		return m


		

	#  
	#          O1
	#       /      \
	#  O2=C3        Ni
	#       \      /
	#        C1==C2
	#       /      \
	#    R1/R2    R2/R1
	#


	#TODO: see if you can find a more elegant solution to grabbing the different ligand indices
	def constructRingIntermediate(self, fiveMemberedRingInstance, switchR1andR2 = False):
		
		fiveMemberedRingInstance.putInPlaneOf_MetalAlkyne(xyzObject = self)
		
		#Add the carboxylate atoms to the xyz string/atom objects
		#NOTE: this is the sanity check
		#Just take the R1 and R2 groups, get their translation vectors first

		if switchR1andR2 == False:
			C1_index, C2_index = self.findX1_andX2nearNickelIndices()
			C1_r, C2_r = self.findX1_andX2nearNickelRs()
		elif switchR1andR2 == True:
			C2_index, C1_index = self.findX1_andX2nearNickelIndices()
			C2_r, C1_r = self.findX1_andX2nearNickelRs()

		else:
			raise Exception("You can only enter True or false for the switchR1andR2 variable")


		#Add a routine for translating C1 and its R group
		if switchR1andR2 == False:
			R1 = self.getRGroup(index = 1, deleteRgroupFromXYZ = False)
			R2 = self.getRGroup(index = 2, deleteRgroupFromXYZ = False)
		
		elif switchR1andR2 == True:
			R1 = self.getRGroup(index = 2, deleteRgroupFromXYZ = False)
			R2 = self.getRGroup(index = 1, deleteRgroupFromXYZ = False)
	
		fmr_C1_r = fiveMemberedRingInstance.C1atom.r
		fmr_C2_r = fiveMemberedRingInstance.C2atom.r

		C1TranslationVector = fmr_C1_r - C1_r
		C2TranslationVector = fmr_C2_r - C2_r		
		
		#Removing the R groups for translating
		self.deleteRgroups()
	
		#Now, do the translation
		R1.translate(C1TranslationVector)
		R2.translate(C2TranslationVector)
		
		#Note, the orientation of the unit vectors produced here are sane
		fmr_C1_H1 = functions.unitVector(fiveMemberedRingInstance.C1atom.r, 
						fiveMemberedRingInstance.R1atom.r)

		fmr_C2_H2 = functions.unitVector(fiveMemberedRingInstance.C2atom.r, 
						fiveMemberedRingInstance.R2atom.r)

		R1.reorientRgroup(fmr_C1_H1)
		R2.reorientRgroup(fmr_C2_H2)

		#I'm grabbing this so my code doesn't need to be agnostic to the scaffold location
		bipyIndices = len(self.atoms)		

		#Now, reinstall the R group to the xyz file
		for atom in R1.atoms:
			self.atoms.append(atom)

		self.R1_indices = list(range(bipyIndices, len(self.atoms)))
		self.R1_Cr = fiveMemberedRingInstance.C1atom.r

		for atom in R2.atoms:
			self.atoms.append(atom)

		self.R2_indices = list(range(self.R1_indices[-1] + 1, len(self.atoms)))
		self.R2_Cr = fiveMemberedRingInstance.C2atom.r

		
		self.atoms.append(fiveMemberedRingInstance.O1atom)
		self.atoms.append(fiveMemberedRingInstance.C3atom)
		self.atoms.append(fiveMemberedRingInstance.O2atom)

		#I decided to add these attributes instead of making my code more elegant
		self.O1atom = fiveMemberedRingInstance.O1atom
		self.C3atom = fiveMemberedRingInstance.C3atom
		self.O2atom = fiveMemberedRingInstance.O2atom
		
		self.regenerateAtomLines(self.atoms)
		
		#So I cann access the scaffold Indices later
		#NOTE: scaffold indices will also contain the R groups
		self.scaffoldIndices = list(range(bipyIndices, len(self.atoms)))

	#NOTE: this code assumes that the five member ring intermediate was generated
	def getCarboxylatedR1andR2s(self, deleteRgroups = True):

		R1_xyzLines = [atom.atomString for i, atom in enumerate(self.atoms) if i in self.R1_indices]
		R2_xyzLines = [atom.atomString for i, atom in enumerate(self.atoms) if i in self.R2_indices]

		

		R1AtomNum = len(R1_xyzLines)
		R2AtomNum = len(R2_xyzLines)

		R1_xyzString = "\n".join(R1_xyzLines)
		R2_xyzString = "\n".join(R2_xyzLines)

		#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#NOTE: I think this is the offending line
		#NOTE: if R1_Cr and R2_Cr were not updated, well, my pooch is screwed
		self.R1group = Rgroup(self.R1_indices, R1_xyzString, self.R1_Cr)
		self.R2group = Rgroup(self.R2_indices, R2_xyzString, self.R1_Cr)

		if deleteRgroups:
			self.deleteRgroups(carboxylated = True)

		return self.R1group, self.R2group

	
	def rotateR(self, R_index = 1, rotationAngle = 5):
		R1, R2 = self.getCarboxylatedR1andR2s()

		if R1.R_xyz == '' or R2.R_xyz == '':
			raise Exception("Looks like your R groups are deleted")

		if R_index == 1:
			R1.rotate(rotationAngle)

		elif R_index == 2:
			R2.rotate(rotationAngle)

		else:
			raise Exception("Your R_index must be 1 or two")


		bipyIndices = len(self.atoms)
		#I'm choosing to simply not update the R groups, a new instance is declared
			#Before I need to worry about them again
		for atom in R1.atoms:
			self.atoms.append(atom)

		for atom in R2.atoms:
			self.atoms.append(atom)

		#I should probably find a better fix than this
		self.atoms.append(self.O1atom)
		self.atoms.append(self.C3atom)
		self.atoms.append(self.O2atom)

		#NOTE: I will most likely need these later
		self.scaffoldIndices = list(range(bipyIndices, len(self.atoms)))
		self.regenerateAtomLines(self.atoms)

	#TODO: flesh out this funciton an you'll be golden
	def reduceRclashes(self, R_index = 1, rotationAngle = 5, debug = False, VMDdebug_type1 = False,
				VMDdebug_type2 = False):

		radianAngle = (math.pi/180)*rotationAngle

		if VMDdebug_type2:
			os.mkdir('./tmp')

		valenceSanityMetrics = []
		Θ = []
		for θ in range(rotationAngle, 360 + rotationAngle, rotationAngle):
			
			self.valenceSanityCheck()

			
			self.rotateR(R_index, radianAngle)

			if VMDdebug_type1:
				self.viewInVMD()
			if VMDdebug_type2:
				self.printToFile(f'./tmp/tmp_{θ}.xyz')
				
			
			if self.valenceSanity == 'sane':
				valenceSanityMetric = 0
			
			elif self.valenceSanity == 'not sane':
				valenceSanityMetric = len(self.valenceSanityRecord.split("\n"))/2

			valenceSanityMetrics.append(valenceSanityMetric)
			Θ.append(θ)

		sanityMinimum = min(valenceSanityMetrics)

		if debug:
			print("........................................")
			print(sanityMinimum)
			print(f"valenceSanityMetrics array {valenceSanityMetrics}")
			print("........................................")


		if VMDdebug_type2:

			subprocess.run(["vmd", "-e", os.path.abspath("loadAllXYZs.tcl")])
			shutil.rmtree('./tmp')
			
		'''
		for valenceSanityMetric, θ in zip(valenceSanityMetrics, Θ):
			if valenceSanityMetric == sanityMinimum:
				optimal_θ = θ + rotationAngle
		'''
		matches = [x == sanityMinimum for x in valenceSanityMetrics]

		optimalIndex = functions.find_middle_of_longest_run(matches)
		optimal_θ = optimalIndex * radianAngle 
		

		self.rotateR(R_index = R_index, rotationAngle = optimal_θ)


	#NOTE: this is just a function for finding the oxygen and carbon nearest
		#NOTE: to the nickel
	'''
	def getC1O1_atoms(self):
			
		self.regenerateAtomLines(self.atoms)
		#Printing to a temp file

		tmpFileName = functions.random_filename(prefix = 'tmp', suffix=".xyz", length=8)
		self.printToFile(tmpFileName)
		atoms = read(tmpFileName)
		#Bad things happen if you don't remove this file immediately after creation
		os.remove(tmpFileName)

		cutoffs = natural_cutoffs(atoms)
		nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
		nl.update(atoms)
	
		#Best add indices
		self.addAtomIndices()
		for atom in self.atoms:
			if atom.element == 'Ni':
				Ni_i = atom.index
				self.metalAtom = atom
				break
	
		indices, _ = nl.get_neighbors(Ni_i)
		for atom in self.atoms:
			if atom.element == 'O' and atom.index in indices:
				self.O1_atom = atom
				break

		for atom in self.atoms:
			if atom.element == 'C' and atom.index in indices:
				self.C1_atom = atom
				break

		return self.C1_atom, self.O1_atom

	
	#NOTE: This defines the nickel attribute, but it is hidden in the code
	def getN1N2_atoms(self):
		
		self.regenerateAtomLines(self.atoms)
		#Printing to a temp file
		tmpFileName = functions.random_filename(prefix = 'tmp', suffix=".xyz", length=8)
		self.printToFile(tmpFileName)
		atoms = read(tmpFileName)
		#Bad things happen if you don't remove this file immediately after creation
		os.remove(tmpFileName)

		cutoffs = natural_cutoffs(atoms)
		nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
		nl.update(atoms)

		#Best add indices
		self.addAtomIndices()
		for atom in self.atoms:
			if atom.element == 'Ni':
				Ni_i = atom.index
				self.metalAtom = atom
				break
	
		indices, _ = nl.get_neighbors(Ni_i)
		for atom in self.atoms:
			if atom.element == 'N' and atom.index in indices:
				self.N1_atom = atom
				break

		for atom in self.atoms:
			if atom.element == 'N' and atom.index in indices and\
				atom.index != self.N1_atom.index:
				self.N2_atom = atom
				break

		return self.N1_atom, self.N2_atom
	'''
	#NOTE: this method was originally designed by yours truly, but I had chatGPT
		#NOTE: design a fix for cases where a nearby hydrogen was registered
		#NOTE: as one of the metal-coordinating atoms
	def getC1O1_atoms(self):
		# make sure strings/indices are in sync
		self.regenerateAtomLines(self.atoms)

		# temp file for ASE
		tmpFileName = functions.random_filename(prefix='tmp', suffix=".xyz", length=8)
		self.printToFile(tmpFileName)
		atoms = read(tmpFileName)
		os.remove(tmpFileName)

		cutoffs = natural_cutoffs(atoms)
		nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
		nl.update(atoms)

		# add indices to our own atoms
		self.addAtomIndices()

		# find Ni
		for atom in self.atoms:
			if atom.element == 'Ni':
				Ni_i = atom.index
				self.metalAtom = atom
				break
		else:
			raise ValueError("No Ni atom found in structure")

		# raw neighbors from ASE (can contain H!)
		neighbor_indices, _ = nl.get_neighbors(Ni_i)

		# turn into real atoms
		neighbor_atoms = [a for a in self.atoms if a.index in neighbor_indices]

		# drop hydrogens that are too close to Ni
		filtered = [a for a in neighbor_atoms if a.element != 'H']

		# sort by distance to Ni so the 'closest real ones' win
		def _dist_to_ni(a):
			return float(np.linalg.norm(a.r - self.metalAtom.r))

		filtered.sort(key=_dist_to_ni)

		# now pick first C and first O from filtered list
		C1_atom = None
		O1_atom = None
		for a in filtered:
			if a.element == 'C' and C1_atom is None:
				C1_atom = a
			elif a.element == 'O' and O1_atom is None:
				O1_atom = a
			if C1_atom is not None and O1_atom is not None:
				break

		if C1_atom is None or O1_atom is None:
			raise RuntimeError(
				f"Could not find both C and O bound to Ni; got neighbors={[(a.element, a.index) for a in filtered]}"
			)

		self.C1_atom = C1_atom
		self.O1_atom = O1_atom
		return self.C1_atom, self.O1_atom


	def getN1N2_atoms(self):
		self.regenerateAtomLines(self.atoms)

		tmpFileName = functions.random_filename(prefix='tmp', suffix=".xyz", length=8)
		self.printToFile(tmpFileName)
		atoms = read(tmpFileName)
		os.remove(tmpFileName)

		cutoffs = natural_cutoffs(atoms)
		nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
		nl.update(atoms)

		self.addAtomIndices()

		# find Ni
		for atom in self.atoms:
			if atom.element == 'Ni':
				Ni_i = atom.index
				self.metalAtom = atom
				break
		else:
			raise ValueError("No Ni atom found in structure")

		neighbor_indices, _ = nl.get_neighbors(Ni_i)
		neighbor_atoms = [a for a in self.atoms if a.index in neighbor_indices]

		# drop H neighbors around the metal
		neighbor_atoms = [a for a in neighbor_atoms if a.element != 'H']

		# sort by distance to metal
		def _dist_to_ni(a):
			return float(np.linalg.norm(a.r - self.metalAtom.r))
		neighbor_atoms.sort(key=_dist_to_ni)

		# pick two N's
		N_atoms = [a for a in neighbor_atoms if a.element == 'N']

		if len(N_atoms) < 2:
			raise RuntimeError(
				f"Expected at least 2 N neighbors to Ni, got {[(a.element, a.index) for a in neighbor_atoms]}"
			)

		self.N1_atom = N_atoms[0]
		self.N2_atom = N_atoms[1]
		return self.N1_atom, self.N2_atom


	def getC1O1_crossProduct(self):
		
		C1atom, O1atom = self.getC1O1_atoms()
		metalAtom = self.metalAtom

		self.C1O1MetalCrossProduct = \
			functions.findCrossProduct(
			C1atom.r, O1atom.r, metalAtom.r)

		return self.C1O1MetalCrossProduct

	def getN1N2_crossProduct(self):
		
		N1atom, N2atom = self.getN1N2_atoms()
		metalAtom = self.metalAtom

		self.N1N2MetalCrossProduct = \
			functions.findCrossProduct(
			N1atom.r, N2atom.r, metalAtom.r)

		return self.N1N2MetalCrossProduct

	#NOTE: this is the dumbest name for a function I've come up with for a long time
		#NOTE: I used vectors for N1 and N2 combined with vectors for C1 and O1
		#NOTE: to get a flavor for how much rotating needs to be done for the five-ring
		#NOTE: intermediate
	def crossProductDebugForTetrahedralAndPlanarSolution(self):
		C1O1_crossProduct = self.getC1O1_crossProduct()
		N1N2_crossProduct = self.getN1N2_crossProduct()
		r = self.metalAtom.r

		C1O1vector = (r, r+C1O1_crossProduct)
		N1N2vector = (r, r+N1N2_crossProduct)

		return C1O1vector, N1N2vector

	def findCrossProductAngles_ofPrimary_SP3(self):

		#Grabbing the cross products
		C1O1_crossProduct = self.getC1O1_crossProduct()
		N1N2_crossProduct = self.getN1N2_crossProduct()	

		#Finding their angle
		θ = functions.findAngleBetweenVectors(C1O1_crossProduct, N1N2_crossProduct)
		
		self.SP3_angle = θ

		return self.SP3_angle

	
	def findMetalAxisOfRotation(self):

		self.findCrossProductAngles_ofPrimary_SP3()
		
		#Collecting N1 and N2
		N1atom, N2atom = self.getN1N2_atoms()
		metalAtom = self.metalAtom

		N1N2_midpoint = functions.findMidpoint(N1atom.r, N2atom.r)

		u_rot = functions.unitVector(N1N2_midpoint, metalAtom.r)

		return u_rot

	#NOTE: Don't forget to regenerate your atom strings when you are done
	def intermediateRotation(self, byAngle = False, angle = 0, debug = False):
		

		u_rot = self.findMetalAxisOfRotation()

		if hasattr(self, "scaffoldIndices"):
			pass
		else:
			raise Exception("You executed the code without constructing your intermediate")

		
		#I think I can avoid copying this data

		scaffoldAtoms = [self.atoms[i] for i in self.scaffoldIndices]
		
		#I'm guessing the code will behave a litle more closely to expectation
			#If I delete some of the atoms then add them later

		self.atoms = [self.atoms[i] for i in range(0,len(self.atoms)) 
						if i not in self.scaffoldIndices]

		P = np.array([scaffoldAtom.r for scaffoldAtom in scaffoldAtoms])
		
		#P_r = functions.rotatePointsByAngle(P, self.metalAtom.r, u_rot, angle)

		if byAngle:
			P_r = functions.rotatePointsByAngle(P, self.metalAtom.r, u_rot, angle)
		else:
			if debug:
				print("/////////////////////////////////////")
				print(f"Number of included atoms {len(P)}")
				print("metal atom positions and cross products of central atoms")
				print(f'self.metalAtom.r: {self.metalAtom.r}')
				print(f'self.metalAtom.r + self.C1O1MetalCrossProduct'
				f'{self.metalAtom.r + self.C1O1MetalCrossProduct}')
				print(f'self.N1N2MetalCrossProduct {self.N1N2MetalCrossProduct}')
				print(f'self.N1N2MetalCrossProduct + metalAtom.r'
					f' {self.N1N2MetalCrossProduct + self.metalAtom.r}')
				print("/////////////////////////////////////")

	
			P_r, antiparallelFlag = functions.reorient_points(P, self.metalAtom.r, self.C1O1MetalCrossProduct, self.N1N2MetalCrossProduct, debug = debug)

			if antiparallelFlag:

				if debug:
					print("The antiparallel code executed")

				P = functions.rotatePointsByAngle(P, self.metalAtom.r, u_rot, angle)
				P_r, _ = functions.reorient_points(P, self.metalAtom.r, self.C1O1MetalCrossProduct, self.N1N2MetalCrossProduct, debug = debug)

				


		for atom, p_r in zip(scaffoldAtoms, P_r):
			atom.r = p_r

		self.atoms+=scaffoldAtoms

		self.regenerateAtomLines(self.atoms)
	
	#Note: seems sane
	def forcePlanar(self, debug = False):
		
		self.intermediateRotation(debug = debug)

	#TODO: needs adjusted to handle radians
	def forceTetrahedral(self):
		#First forcing planar:
		self.intermediateRotation()

		angleInRadians = (90/180)*math.pi

		#Then forcing 90°
		self.intermediateRotation(byAngle = True, angle = angleInRadians)

	def pivotIntermediate(self, angleInDegrees = 5):

		angleInRadians = (angleInDegrees/180)*math.pi

		self.intermediateRotation(byAngle = True, angle = angleInRadians)

	#I do not mind the weight because I see no point
	def findCOMproxy(self):
		
		R = [atom.r for atom in self.atoms]
		R = np.array(R)

		COM = functions.findCentralPoint(R)
		
		self.COM = COM
		return self.COM

	#TODO: draft this method
	#NOTE: Implement this after generating the five member ring scaffold
	def separateFiveMemberRing(self, fiveMemberedRingInstance, d_sep):


		bipy_COM = self.findCOMproxy()
		fmr_COM = fiveMemberedRingInstance.findCOMproxy()

		v_fmr_bpy = functions.unitVector(fmr_COM, bipy_COM) * d_sep
		v_bpy_fmr = functions.unitVector(bipy_COM, fmr_COM) * d_sep

		try:
			bipyIndices = list(range(0, self.scaffoldIndices[0]))
			scaffoldIndices = list(range(self.scaffoldIndices[0], self.scaffoldIndices[-1]))

		except:
			raise Exception("Most probably, you have neglected to execute the .constructRingIntermediate method yet")

		for i, atom in enumerate(self.atoms):

			if i in bipyIndices:
				atom.r+=v_fmr_bpy
				atom.updateString()

			if i in scaffoldIndices:
				atom.r+=v_bpy_fmr
				atom.updateString()
		
		self.regenerateAtomLines(self.atoms)
	
	def addAtomIndices(self):
		#Regenerating atoms does not harm
		self.regenerateAtomLines(self.atoms)

		for i, atom in enumerate(self.atoms):
			atom.index = i



	def addAtomNeighborNumberToAtoms(self):
		
		#Regenerating atoms does no harm
		self.regenerateAtomLines(self.atoms)
		#Atom indices will be included by default
		self.addAtomIndices()
		#Printing to a temp file
		tmpFileName = functions.random_filename(prefix = 'tmp', suffix=".xyz", length=8)
		self.printToFile(tmpFileName)
		#Running the ASE environment
		cutoffs = natural_cutoffs(read(tmpFileName))
		nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
		nl.update(read(tmpFileName))
		#Removing the file or shenanigans begin
		os.remove(tmpFileName)

		for atom in self.atoms:	
			indices, _ = nl.get_neighbors(atom.index)
			atom.neighborIndices = indices

	#TODO: battle test this function
	def valenceSanityCheck(self):
		#You should access your variables with the MAX_VALENCE dictionary

		self.addAtomNeighborNumberToAtoms()

		#Innocent untilProven guilty
		self.valenceSanity = 'sane'

	
		for atom in self.atoms:
			if MAX_VALENCE[atom.element] < len(atom.neighborIndices):
				self.valenceSanity = 'not sane'
				atom.valenceSanity = 'not sane'
			else:
				atom.valenceSanity = 'sane'

		self.valenceSanityRecord = ''

		if self.valenceSanity == 'not sane':
			for atom in self.atoms:
				if atom.valenceSanity == 'not sane':
					self.valenceSanityRecord += f"Atom of element {atom.element}\n"\
							f"@ index {atom.index} has abnormal valence \n"

		else:
			self.valenceSanityRecord = ''
			self.valenceSanity = 'not sane'

	#TODO: the structures being loaded into my sanity check don't seem to be updated

	
	def pivotCorrectionForValenceSanity(self, rotationAngle=5, debug = False):

		#This one already accepts angles in degrees
		rotationCount = 0

		# Ensure strings are in sync before we start
		self.regenerateAtomLines(self.atoms)

		valenceSanityMetrics = []
		Θ = []
		

		for θ in range(rotationAngle, 360 + rotationAngle, rotationAngle):

			# 1) pivot first (updates atom.r; intermediateRotation ends with regenerate)
			self.pivotIntermediate(angleInDegrees=rotationAngle)
			
			# 2) I need to keep track of rotations in practice
			rotationCount+=1

			# 3) regenerate strings (belt-and-suspenders)
			self.regenerateAtomLines(self.atoms)

			# 4) re-check using the updated geometry
			self.valenceSanityCheck()
			
			# 5) quantifying sanity
			if self.valenceSanityRecord == '':
				valenceSanityMetric = 0
			
			else:
				valenceSanityMetric = len(self.valenceSanityRecord.split("\n"))/2
				
			#6) Adding the metrics right where I need it
			valenceSanityMetrics.append(valenceSanityMetric)
			Θ.append(θ)
			
			if debug == True:
				print("********************************************")
				print(f'{rotationCount*rotationAngle} degrees')  # rotations applied
				print(self.valenceSanityRecord)  # this now reflects post-pivot bonds
				self.viewInAvogadro()
				# OPTIONAL: visualize current geometry (post-pivot)
				print("********************************************")

		
		sanityMinimum = min(valenceSanityMetrics)

		for i, valenceSanityMetric in enumerate(valenceSanityMetrics):
			if valenceSanityMetric == sanityMinimum:
				optimal_θ = rotationAngle*(i + 1)
				break
		
		self.pivotIntermediate(angleInDegrees=optimal_θ)
			

	def writeSanityRecord(self, sanityRecordPath):
		
		self.valenceSanityCheck()

		if self.valenceSanityRecord == '':
			return 0
		else:
			file = open(sanityRecordPath, "a")
			file.write(self.valenceSanityRecord)
			file.close()
	
#NOTE: R_xyz are the lines, not the complete strings to be seen in xyz files 
class Rgroup:
	def __init__(self, R_indices, R_xyz, C_r):
		self.R_indices = R_indices
		self.R_xyz = R_xyz
		#Note: this attribute contains some lines from the xyz file
			#Note: it does not contain all the information needed for a .xyz
		self.R_xyzFileContents = str(len(self.R_indices)) + "\n" + "\n" + self.R_xyz
		self.atoms = [atom(line) for line in self.R_xyz.split("\n")]
		self.atomLines = [line for line in self.R_xyz.split("\n")]
		
		self.C_r = C_r.copy()
		self.origin = C_r.copy()

	#NOTE: updatge for both the atoms and R_xyz attribute
	def translate(self, vectorOfTranslation):
		self.R_xyz = ''
		for atom in self.atoms:
			atom.translate(vectorOfTranslation)
			self.R_xyz += f"{atom.atomString}\n"
		self.C_r += vectorOfTranslation
		self.origin += vectorOfTranslation

	def getRgroupOrientation(self):

		#displacements = [atom.r - self.origin for atom in self.atoms]
		displacements = []
		for atom in self.atoms:
			displacements.append(atom.r - self.origin)
		avgVector = np.mean(displacements, axis=0)
		norm = np.linalg.norm(avgVector)
		return avgVector / norm


	#NOTE: the update takes place within the R object
	def reorientRgroup(self, newOrientation):
		#Note: this returns non-type figure out why
		R_u = self.getRgroupOrientation()

		R_P = np.stack([atom.r for atom in self.atoms])
		

		#NOTE: I am worried these points are not getting updated
		points_1, _ = functions.reorient_points(R_P, self.origin, R_u, newOrientation)
		for atom, point in zip(self.atoms, points_1):
			atom.r = point
			atom.updateString()

	def regenerateAtomLines(self, atoms):
		
		for atom in self.atoms:
			atom.updateString()

		self.atomCount = len(atoms)
		self.atomLines = [atom.atomString for atom in atoms]
		self.R_xyz = "\n".join(self.atomLines)
		self.R_xyzFileContents = str(self.atomCount) + "\n" + "\n" + "\n".join(self.atomLines)


	def printToFile(self, xyzPath):
		file = open(xyzPath, "a")
		file.write(self.R_xyzFileContents)
		file.close()

	def viewInVMD(self):
		# Courtesy of chatGPT
        	with tempfile.NamedTemporaryFile(delete=False, suffix=".xyz") as tmp:
            		tmp.write(self.R_xyzFileContents.encode("utf-8"))
            		tmp_path = tmp.name

        	try:
            		subprocess.run(["vmd", tmp_path])

        	finally:
            		if os.path.exists(tmp_path):
                		os.remove(tmp_path)

	def updateOrigin(self):

		self.C_r = self.atoms[0].r
		self.origin = self.atoms[0].r

	def rotate(self, angleInDegrees = 5):
		self.updateOrigin()

		u_rot = self.getRgroupOrientation()
		
		P = np.array([atom.r for atom in self.atoms])

		P_r = functions.rotatePointsByAngle(P, self.origin, u_rot, angleInDegrees)

		for atom, p_r in zip(self.atoms, P_r):
			atom.r = p_r

		self.regenerateAtomLines(self.atoms)

			
	

			
		
		
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

#  
#          O1
#       /      \
#  O2=C3        Ni
#       \      /
#        C1==C2
#       /      \
#    R1/R2    R2/R1
#

        
class fiveMemberedRing:
	def __init__(self, pathToXYZfile = "/home/alpal/projects/methanCapture/carboxylationProblem/data/carboxylatedExample.xyz"):

		xyzStructureString = open(pathToXYZfile).read()
		xyzStructureInstance = xyzStructure(xyzStructureString)

                        
		self.O1Index = 26
		self.C3Index = 21
		self.O2Index = 27
		self.C1Index = 22
		self.C2Index = 20
		self.NiIndex = 25
		self.R1Index = 24
		self.R2Index = 23
        
		self.d_O1C3 = xyzStructureInstance.findAtom_AtomDistance(self.O1Index, self.C3Index)
		self.d_C3O2 = xyzStructureInstance.findAtom_AtomDistance(self.C3Index, self.O2Index)
		self.d_C3C1 = xyzStructureInstance.findAtom_AtomDistance(self.C3Index, self.C1Index)
		self.d_C1C2 = xyzStructureInstance.findAtom_AtomDistance(self.C1Index, self.C2Index)
		self.d_C2Ni = xyzStructureInstance.findAtom_AtomDistance(self.C2Index, self.NiIndex)
		self.d_NiO1 = xyzStructureInstance.findAtom_AtomDistance(self.NiIndex, self.O1Index)

		self.O1atom = xyzStructureInstance.returnAtomByIndex(self.O1Index)
		self.C3atom = xyzStructureInstance.returnAtomByIndex(self.C3Index)
		self.O2atom = xyzStructureInstance.returnAtomByIndex(self.O2Index)
		self.C1atom = xyzStructureInstance.returnAtomByIndex(self.C1Index)
		self.C2atom = xyzStructureInstance.returnAtomByIndex(self.C2Index)
		self.Niatom = xyzStructureInstance.returnAtomByIndex(self.NiIndex)
		self.R1atom = xyzStructureInstance.returnAtomByIndex(self.R1Index)
		self.R2atom = xyzStructureInstance.returnAtomByIndex(self.R2Index)

		self.atoms = [self.O1atom, self.C3atom, self.O2atom, self.C1atom, 
				self.C2atom, self.Niatom, self.R1atom, self.R2atom]



	#NOTE: this function is untested
	#NOTE: this function will give me the points, but it will not place atoms in my xyzObject
		#NOTE: instance
	#NOTE: I think this should enforce sane geometries on its own
		#NOTE: I can only verify after adding
	

	def putInPlaneOf_MetalAlkyne(self, xyzObject):

		src_anchor = np.vstack([self.O1atom.r, self.Niatom.r, self.C2atom.r])

		C1_r, C2_r = xyzObject.findX1_andX2nearNickelRs()
		Ni_r = xyzObject.findNickelR()
		dst_anchor = np.vstack([np.array(C1_r), np.array(Ni_r), np.array(C2_r)])

		R, t = functions.rigid_transform(src_anchor, dst_anchor)

		ring_atoms = [
			self.O1atom,
			self.C3atom,
			self.O2atom,
			self.C1atom,
			self.C2atom,
			self.Niatom,
			self.R1atom,
			self.R2atom
		]

		for atom in ring_atoms:
			atom.r = R @ atom.r + t
			atom.updateString()

	def findCOMproxy(self):
		
		R = [atom.r for atom in self.atoms]
		R = np.array(R)

		COM = functions.findCentralPoint(R)
		
		self.COM = COM
		return self.COM

#TODO: add a method to the xyzStructure for identifying the indices of those atoms which are attached to C1 and C2 (this could be a headache..... I'll probably need a molecular graphics software to guess bonds for me

#TODO: add an R1andR2 object for handling rotation operations and such
	#TODO: it is best if you can add these as attributes to the xyzStructure

#TODO: add a method for paritionioning R groups from the xyz file

#TODO: add a method for fishing out bond lengths from the example structure you constructed in gaussview

#TODO: add a method for putting the five member ring in the plane of the metal alkyne
#TODO: add a method for distancing C1 and C2 from one another
#TODO: add a method to the carboxylate method for fixing bond lengths
#TODO: add a method to the carboxylate for adding a hydrogen when you configure it to handle detachment
	#TODO: from nickel
#TODO: add a method for finding the angles between two atoms to handle rotations

#TODO: for whatever reason, the fiveMember ring attributes are not
			#TODO: refelecting the transformation

#TODO: battle test this code by adding the five member ring

#TODO: rename the atoms attribute, it makes your code confusing

#TODO: Build a force-field optimization method (i.e. UFF or MMFF), make sure it doesn't distort the
	#TODO: geometry significantly

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#TODO: Try an RDkit-based method for bond-length sanity checks
	#NOTE: I'm not so sure this is reasonable because:
	#NOTE: 1 ASE assigns bonds based on presumably sane criterion, it's not like
	#NOTE: a topolgo file where you can hypothetically force an insane bond length by mistake
	#NOTE: 2: DFT or DFT analogs (i.e. machine learning solutions) will most probably adjust
	#NOTE: bond lengths
	#NOTE: the priority seems to be valence sanity checks
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#TODO: finish this function

#TODO: finish fleshing out this function

	
'''
	#NOTE: Looks sane
	def generateRDKitMol(self):
		self.printToFile('tmp.xyz')
		mol = MolFromXYZFile('tmp.xyz')
		os.remove('tmp.xyz')
		return mol

	def debug(self):
		raw_mol = self.generateRDKitMol()
		mol = Chem.Mol(raw_mol)
		rdDetermineBonds.DetermineBonds(mol)
		raise Exception("Stop here")
		#for atom in mol.GetAtoms():
		#	print(atom.GetDegree())
		
		#rdmolops.Kekulize(mol, clearAromaticFlags=True)
		
		for atom in mol.GetAtoms():
			print(atom.GetDegree())

		#mol.Debug()

'''

