import tempfile
import subprocess
import os

class xyzStructure:
	def __init__(self, xyzString):
		self.xyzString = xyzString

        def findC1andC2

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

	def findX1_andX2nearNickelIndices(self, prefer_two_shortest=True, scale=1.25, metal_fudge=0.20, particleType="c"):
		xzy_text = self.putInXYZfileFormat()
		atoms = parse_xyz(xzy_text)

		if not atoms:
			raise ValueError("No atoms parsed from XYZ.")

		symbols = [s for s, _ in atoms]
		positions = [p for _, p in atoms]
		adj = build_connectivity(symbols, positions, scale=scale, metal_fudge=metal_fudge)

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
			c_neighbors_sorted = sorted(c_neighbors, key=lambda j: dist(positions[i_ni], positions[j]))

			# NOTE: c_neighbors_sorted does not look empty
			if prefer_two_shortest:
				C1, C2 = c_neighbors_sorted[:2]
			else:
				C1, C2 = c_neighbors_sorted[0], c_neighbors_sorted[1]

			# Ensure consistent ordering by index (optional)
			if C2 < C1:
				C1, C2 = C2, C1

			return C1, C2

	def findX1_andX2nearNickelRs(self, particleType="c"):
		C1_i, C2_i = self.findX1_andX2nearNickelIndices(particleType=particleType)

		C1 = self.muratsXyzTuple[C1_i]
		C2 = self.muratsXyzTuple[C2_i]

		C1_r = C1[-1]
		C2_r = C2[-1]

		return C1_r, C2_r


