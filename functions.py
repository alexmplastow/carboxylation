import sqlite3
import pandas as pd
import json
import math
from collections import defaultdict
import numpy as np
import objects

def getPandasDFfromDB(pathToDBfile):
	#My friend Murat wrote these two lines, I have no idea what the squlite3 routine is doing
	conn = sqlite3.connect(pathToDBfile)
	database_df = pd.read_sql("SELECT * FROM calculations", conn)
	return database_df

def getXYZstructureList(pathToDBfile):
	df = getPandasDFfromDB(pathToDBfile)
	optimized_xyz_list = []
	for index, row in df.iterrows():
		try:
			blob_data_dict = json.loads(row['blob_data'])
			optimized_xyz_list.append(blob_data_dict['opt_xyz'])
		except (json.JSONDecodeError, KeyError) as  e:
			print(f"Error processing row {index}: {e}")
			optimized_xyz_list.append(None)
	return optimized_xyz_list


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

def bonded(sym_i, sym_j, d, scale=1.25, metal_fudge=0.20):
	"""
	Heuristic bond criterion:
	  d <= scale * (r_i + r_j) + (metal_fudge if either is a metal like Ni)
	"""
	base = scale * (rcov(sym_i) + rcov(sym_j))
	extra = metal_fudge if (sym_i.lower() in ("ni",) or sym_j.lower() in ("ni",)) else 0.0
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

#Courtesy of chatGPT
def rotatePoints(points, origin, u, v):
    """
    Rotate a set of 3D points so that orientation vector u aligns with v.

    Parameters
    ----------
    points : (N, 3) array
        3D coordinates of points
    origin : (3,) array
        The origin about which to rotate
    u, v : (3,) arrays
        Unit vectors: current orientation and target orientation

    Returns
    -------
    rotated_points : (N, 3) array
    """
    u = u / np.linalg.norm(u)
    v = v / np.linalg.norm(v)

    # Check if vectors are parallel
    if np.allclose(u, v):
        return points.copy()
    if np.allclose(u, -v):
        # Rotate 180° around any axis perpendicular to u
        perp = np.array([1,0,0]) if abs(u[0]) < 0.9 else np.array([0,1,0])
        k = np.cross(u, perp)
        k /= np.linalg.norm(k)
        theta = np.pi
    else:
        k = np.cross(u, v)
        k /= np.linalg.norm(k)
        theta = np.arccos(np.dot(u, v))

    # Rodrigues rotation formula
    K = np.array([
        [0, -k[2], k[1]],
        [k[2], 0, -k[0]],
        [-k[1], k[0], 0]
    ])
    R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)

    # Apply rotation
    shifted = points - origin
    rotated = shifted @ R.T
    return rotated + origin

	


