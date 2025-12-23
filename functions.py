import sqlite3
import pandas as pd
import json
import math
from collections import defaultdict
import numpy as np
import objects
import random
import string
import re
import sqlite3
import ast

def getPandasDFfromDB_jsonFix(pathToDBfile):

    conn = sqlite3.connect(pathToDBfile)
    df = pd.read_sql("SELECT * FROM calculations", conn)
    conn.close()

    # Convert JSON strings in 'blob_data' column to dicts
    df['blob_data'] = df['blob_data'].apply(json.loads)
    return df

'''
def getPandasDFfromDB(pathToDBfile):
	#My friend Murat wrote these two lines, I have no idea what the squlite3 routine is doing
	conn = sqlite3.connect(pathToDBfile)
	database_df = pd.read_sql("SELECT * FROM calculations", conn)
	return database_df
'''

def getPandasDFfromDB(pathToDBfile):
    conn = sqlite3.connect(pathToDBfile)
    database_df = pd.read_sql("SELECT * FROM calculations", conn)

    # Convert 'blob_data' strings to dictionaries
    if 'blob_data' in database_df.columns:
        def parse_blob(s):
            if not isinstance(s, str):
                return s
            try:
                return json.loads(s)
            except json.JSONDecodeError:
                try:
                    return ast.literal_eval(s)
                except Exception:
                    return s  # fallback: leave as-is if parsing fails

        database_df['blob_data'] = database_df['blob_data'].apply(parse_blob)

    return database_df

def getXYZstructureList(pathToDBfile):

	df = getPandasDFfromDB(pathToDBfile)
	optimized_xyz_list = []

	for index, row in df.iterrows():
		optimized_xyz_list.append(row['blob_data']['opt_xyz'])
	return optimized_xyz_list

def load_xyz_structures_from_file(filepath):
    """
    Reads a text file containing multiple XYZ structures concatenated together,
    and returns a list of xyzStructure instances.
    """
    with open(filepath, "r") as f:
        lines = [ln.rstrip() for ln in f if ln.strip()]

    structures = []
    i = 0
    while i < len(lines):
        try:
            atom_count = int(lines[i])
        except ValueError:
            i += 1
            continue  # skip non-count lines (like stray atoms or headers)

        # Build the XYZ block: count + name + atom lines
        block_end = i + 2 + atom_count
        if block_end > len(lines):
            break  # incomplete block at end of file

        xyz_block = "\n".join(lines[i:block_end])
        structures.append(objects.xyzStructure(xyz_block))
        i = block_end  # move to next structure

    return structures

def getXYZandBuildingBlock_structureList(pathToDBfile, xyzFile, sanityCheck = False):
	
	df = getPandasDFfromDB(pathToDBfile)
	xyzIndices = []
	for i in df['blob_data']:
		#I think this value is the corresponding structure
		xyzIndex = int(i['unique_name'].split("_")[3])
		xyzIndices.append(xyzIndex)

	xyzInstances = load_xyz_structures_from_file(xyzFile)
		

	structureCodes = [xyzInstances[i_x].name_line for i_x in xyzIndices]

	optimized_xyz_list = []
	initial_xyz_list = []


	for index, row in df.iterrows():
		initial_xyz_list.append(row['blob_data']['initial_xyz'])
		optimized_xyz_list.append(row['blob_data']['opt_xyz'])

	if sanityCheck:
		commonIndex = 0
		#Checking the first XYZ instance from optimized XYZ
		print("****************************")
		print(initial_xyz_list[commonIndex])
		print("///////////////////////////////")
		print(xyzInstances[xyzIndices[commonIndex]].xyzString)
		print("///////////////////////////////")
		print(structureCodes[commonIndex])
		print("****************************")
		raise Exception("Stop here")


	return optimized_xyz_list, structureCodes



#Written by Murat (but probably by chatGPT or Gemini)
def parse_xyz(xyz_text):
	"""
	Parse an XYZ string into a list of (symbol, (x,y,z)).
	Works whether the file has:
	  - no header
	  - atom count only
	  - atom count + (optional) comment line (possibly empty)
	"""
	raw_lines = xyz_text.splitlines()
	# keep empty lines to detect presence/absence of comment line properly
	lines = [l.rstrip("\n") for l in raw_lines]

	# Trim leading/trailing blank lines for robustness
	while lines and not lines[0].strip():
		lines.pop(0)
	while lines and not lines[-1].strip():
		lines.pop()

	start = 0
	if lines:
		first = lines[0].strip()
		try:
			nat = int(first)
			# If there is a second line and it *doesn't* look like an atom line,
			# treat it as a comment and skip it; otherwise, assume no comment.
			if len(lines) >= 2:
				parts = lines[1].split()
				if len(parts) < 4 or not parts[0].isalpha():
					start = 2  # count + comment
				else:
					start = 1  # count only (no comment line)
			else:
				start = 1
		except ValueError:
			start = 0

	atoms = []
	for line in lines[start:]:
		if not line.strip():
			continue
		parts = line.split()
		if len(parts) < 4:
			continue
		sym = parts[0]
		try:
			x, y, z = map(float, parts[1:4])
		except ValueError:
			continue
		atoms.append((sym, (x, y, z)))
	return atoms

def dist(a, b):
    return math.sqrt(sum((ai - bi) ** 2 for ai, bi in zip(a, b)))

'''
def bonded(sym_i, sym_j, d, scale=1.25, metal_fudge=0.20):
	"""
	Heuristic bond criterion:
	  d <= scale * (r_i + r_j) + (metal_fudge if either is a metal like Ni)
	"""
	base = scale * (rcov(sym_i) + rcov(sym_j))
	extra = metal_fudge if (sym_i.lower() in ("ni",) or sym_j.lower() in ("ni",)) else 0.0
	return d <= (base + extra)

'''

def bonded(sym_i, sym_j, d, scale=1.25, metal_fudge=0.20):
    """
    Heuristic bond criterion:
      d <= scale * (r_i + r_j) + (metal_fudge if either is a metal like Ni)
    Prevents H from bonding to metal centers.
    """
    # Ignore spurious Ni–H or other metal–H "bonds"
    metals = {"ni", "fe", "cu", "co", "zn", "pd", "pt", "ru", "rh"}
    if (sym_i.lower() in metals and sym_j.lower() == "h") or \
       (sym_j.lower() in metals and sym_i.lower() == "h"):
        return False

    base = scale * (rcov(sym_i) + rcov(sym_j))
    extra = metal_fudge if (sym_i.lower() in metals or sym_j.lower() in metals) else 0.0
    return d <= (base + extra)


COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
    "Ni": 1.24,  # approximate
    # add more if needed
}

def rcov(sym):
	return COVALENT_RADII.get(sym.capitalize(), 0.76)

def build_connectivity(symbols, positions, scale=1.25, metal_fudge=0.20):
	"""Return adjacency list inferred from distances."""
	n = len(symbols)
	adj = defaultdict(set)
	for i in range(n):
		for j in range(i+1, n):
			d = dist(positions[i], positions[j])
			if bonded(symbols[i], symbols[j], d, scale=scale, metal_fudge=metal_fudge):
				adj[i].add(j)
				adj[j].add(i)
	return adj


#Courtesy of chatGPT
def rigid_transform(src_pts, dst_pts):
	"""
	Compute rigid transform that aligns src_pts to dst_pts.
	Each is a (3,3) array with rows as points.
	"""
	src_centroid = np.mean(src_pts, axis=0)
	dst_centroid = np.mean(dst_pts, axis=0)

	# Center
	src_centered = src_pts - src_centroid
	dst_centered = dst_pts - dst_centroid

	# Rotation via SVD (Kabsch algorithm)
	H = src_centered.T @ dst_centered
	U, S, Vt = np.linalg.svd(H)
	R = Vt.T @ U.T
	if np.linalg.det(R) < 0:
		Vt[-1, :] *= -1
		R = Vt.T @ U.T

	# Translation
	t = dst_centroid - R @ src_centroid
	return R, t

def findCrossProduct(q1, q2, q3):
	v1 = q3 - q1
	v2 = q2 - q1
	cross = np.cross(v1, v2)
	return cross/np.linalg.norm(cross)

#Note: based on rodrigues' rotation formula
def rotatePointsAboutCrossProduct(points, cross, angle, pivot_index=0):

    pivot = points[pivot_index]
    shifted = points - pivot

    K = np.array([
        [0, -cross[2], cross[1]],
        [cross[2], 0, -cross[0]],
        [-cross[1], cross[0], 0]
    ])

    R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)

    rotated_shifted = shifted @ R.T
    rotated = rotated_shifted + pivot

    return rotated

def unitVector(point1, point2):
    v = np.array(point2) - np.array(point1)
    return v / np.linalg.norm(v)

#Courtesy of chatGPT
def orientVector(points):
    """
    Given an array of 3D points (shape: N x 3), return
    the dominant orientation as a unit vector.
    """
    # Subtract centroid to center the data
    centroid = np.mean(points, axis=0)
    centered = points - centroid

    # Covariance matrix
    cov = np.cov(centered.T)

    # Eigen decomposition
    eigvals, eigvecs = np.linalg.eig(cov)

    # Pick eigenvector with largest eigenvalue
    principal_axis = eigvecs[:, np.argmax(eigvals)]

    # Normalize to unit vector
    unit_vec = principal_axis / np.linalg.norm(principal_axis)

    return unit_vec

#Takes angles in radians, it sucks
def rotation_matrix(axis, angle):
    """
    Rodrigues' rotation formula for 3D rotation matrix.
    """
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    I = np.eye(3)
    return I + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)

def rotate_points(points, origin, axis, angle_deg):
    """
    Rotate a set of 3D points around an axis through a given origin.
    
    points: (N,3) array of 3D points
    origin: (3,) array, the pivot point
    axis: (3,) array, axis of rotation
    angle_deg: float, rotation angle in degrees
    """
    R = rotation_matrix(axis, angle_deg)
    points = np.asarray(points, dtype=float)
    origin = np.asarray(origin, dtype=float)

    # Translate points so origin is at (0,0,0)
    shifted = points - origin
    # Apply rotation
    rotated = shifted @ R.T
    # Translate back
    return rotated + origin

#NOTE: I think this function would be rather handy
'''
def reorient_points(points, origin, current_dir, target_dir):
    """
    Rotate a set of points so that `current_dir` aligns with `target_dir`.

    Parameters
    ----------
    points : (N,3) array-like
        The set of 3D points to transform.
    origin : (3,) array-like
        The origin point (remains fixed during transformation).
    current_dir : (3,) array-like
        Unit vector giving current orientation of the point set.
    target_dir : (3,) array-like
        Unit vector giving desired orientation of the point set.

    Returns
    -------
    transformed_points : (N,3) ndarray
        The reoriented set of points.
    """

    points = np.array(points, dtype=float)
    origin = np.array(origin, dtype=float)
    u = np.array(current_dir, dtype=float)
    v = np.array(target_dir, dtype=float)

    # Normalize to be safe
    u /= np.linalg.norm(u)
    v /= np.linalg.norm(v)

    # Compute rotation axis (cross product) and angle
    axis = np.cross(u, v)
    axis_norm = np.linalg.norm(axis)

    if axis_norm < 1e-10:
        # u and v are parallel (or antiparallel)
        if np.dot(u, v) > 0:
            R = np.eye(3)  # no rotation needed
        else:
            # 180° rotation: choose any orthogonal axis
            axis = np.array([1.0, 0.0, 0.0])
            if np.abs(u[0]) > 0.9:  # avoid collinearity
                axis = np.array([0.0, 1.0, 0.0])
            axis = axis - np.dot(axis, u) * u
            axis /= np.linalg.norm(axis)
            angle = np.pi
            R = rotation_matrix(axis, angle)
    else:
        axis /= axis_norm
        angle = np.arccos(np.clip(np.dot(u, v), -1.0, 1.0))
        R = rotation_matrix(axis, angle)

    # Apply rotation about the origin
    shifted_points = points - origin
    rotated_points = shifted_points @ R.T
    transformed_points = rotated_points + origin

    return transformed_points
'''

def rotation_matrix(axis, angle):
    """
    Rodrigues' rotation formula for a 3D rotation matrix.
    Ensures the axis is normalized for numerical stability.
    """
    axis = np.asarray(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    K = np.array([
        [0, -axis[2], axis[1]],
        [axis[2], 0, -axis[0]],
        [-axis[1], axis[0], 0]
    ])
    I = np.eye(3)
    return I + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)

def reorient_points(points, origin, current_dir, target_dir, debug=False, antiparallel_threshold=-0.90):
    """
    Rotate a set of 3D points so that `current_dir` aligns with `target_dir`.

    Returns
    -------
    transformed_points : (N,3) ndarray
        The rotated coordinates.
    antiparallel_flag : bool
        True if the vectors were nearly antiparallel (dot < threshold).
    """
    points = np.asarray(points, dtype=float)
    origin = np.asarray(origin, dtype=float)
    u = np.asarray(current_dir, dtype=float)
    v = np.asarray(target_dir, dtype=float)

    # === Normalize and check validity ===
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu < 1e-8 or nv < 1e-8:
        raise ValueError(f"Zero-length direction vector(s): |u|={nu}, |v|={nv}")
    u /= nu
    v /= nv

    dot = np.clip(np.dot(u, v), -1.0, 1.0)
    cross = np.cross(u, v)
    norm_cross = np.linalg.norm(cross)
    angle = math.atan2(norm_cross, dot)

    # === Antiparallel detection ===
    antiparallel_flag = dot < antiparallel_threshold

    if debug:
        print("\n--- reorient_points diagnostics ---")
        print(f"current_dir (u): {u}")
        print(f"target_dir  (v): {v}")
        print(f"dot(u,v): {dot:.6f}")
        print(f"angle (deg): {math.degrees(angle):.6f}")
        print(f"antiparallel_flag: {antiparallel_flag}")
        print("------------------------------------")

    # === Rotation matrix selection ===
    if dot > 0.999999:
        R = np.eye(3)  # nearly parallel
    elif antiparallel_flag:
        # placeholder stable rotation (π about orthogonal axis)
        axis = np.cross(u, [1, 0, 0])
        if np.linalg.norm(axis) < 1e-3:
            axis = np.cross(u, [0, 1, 0])
        axis /= np.linalg.norm(axis)
        R = rotation_matrix(axis, np.pi)
    else:
        axis = cross / norm_cross
        R = rotation_matrix(axis, angle)

    # === Apply rotation ===
    shifted = points - origin
    rotated = shifted @ R.T
    transformed = rotated + origin

    return transformed, antiparallel_flag

def findCentralPoint(points):
	
	N = len(points)
	center = np.array([0, 0, 0], dtype = 'float64')

	for r in points:
		center+=r

	COM = center/N
	return COM

#Courtesy of chatGPT
def findAngleBetweenVectors(v1, v2):
	# Compute dot product and magnitudes
	dot_product = np.dot(v1, v2)
	norm_product = np.linalg.norm(v1) * np.linalg.norm(v2)

	if norm_product == 0:
		raise ValueError("One of the vectors has zero magnitude.")

	# Clip value to avoid numerical issues outside [-1,1]
	cos_theta = np.clip(dot_product / norm_product, -1.0, 1.0)
	angle = np.arccos(cos_theta)

	return np.degrees(angle)
	
def findMidpoint(r1, r2):
	m = (r1 + r2) / 2
	return m

#NOTE: utilizes the rotation matrix function 
#Courtesy of chatGPT
def rotatePointsByAngle(points, origin, axis, angle_deg):
    """
    Rotate a set of 3D points around an axis through a given origin.
    
    points: (N,3) array of 3D points
    origin: (3,) array, the pivot point
    axis: (3,) array, axis of rotation
    angle_deg: float, rotation angle in degrees
    """
    R = rotation_matrix(axis, angle_deg)

    # Translate points so origin is at (0,0,0)
    shifted = points - origin
    # Apply rotation
    rotated = shifted @ R.T
    # Translate back
    return rotated + origin

def find_middle_of_longest_run(matches):
	"""
	Given a list of 0s and 1s, find the index of the middle element
	of the longest continuous run of 1s.

	If there are multiple equally long runs, the first is chosen.
	If no 1s are present, returns None.
	"""
	best_start, best_len = -1, 0
	current_start, current_len = -1, 0

	for i, val in enumerate(matches):
		if val == 1:
			if current_len == 0:  # new run starts
				current_start = i
			current_len += 1
			# update best run if current is longer
			if current_len > best_len:
				best_start, best_len = current_start, current_len
		else:
			current_len = 0  # reset on 0

	if best_len == 0:
		return None  # no runs at all

	middle_index = best_start + best_len // 2
	return middle_index

#Courtesy of chatGPT
def random_filename(prefix='tmp', suffix=".xyz", length=8):
	"""
	Generate a random temporary filename with the given prefix and suffix.
	Does not create the file; only returns a unique name.
	"""
	random_part = ''.join(random.choices(string.ascii_lowercase + string.digits, k=length))
	return f"{prefix}_{random_part}{suffix}"

def extract_ligand_numbers(file_list):
    """
    Extract numeric indices from filenames like 'ligands_6180.xyz'.
    Works even if some filenames omit the 'ligands_' prefix.
    """
    numbers = []
    for f in file_list:
        match = re.search(r"(\d+)\.xyz$", f)
        if match:
            numbers.append(int(match.group(1)))
    return numbers



