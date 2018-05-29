#!/usr/bin/env python

"""
Routines to calculate the Accessible Surface Area of a set of atoms.
The algorithm is adapted from the Rose lab's chasa.py, which uses
the dot density technique found in:

Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent
of Protein Atoms. Lysozyme and Insulin." JMB (1973) 79:351-371.
"""

import math
from vector3d import pos_distance, Vector3d, pos_distance_sq
import numpy as np



#@profile
def find_neighbor_indices(points, radii, probe, k):
    """
    Returns list of indices of atoms within probe distance to atom k. 
    """
    radius = radii[k] + probe + probe
    test_radii = (radius + radii) ** 2

    dist_sq = np.sum(((points - points[k,:])** 2), axis=1)
    
    neighbor_indices = (dist_sq < test_radii)
 
    #Must remove self distance
    neighbor_indices[k] = False

    #Need to return a list/array of integers not a boolean array
    #return (np.where(neighbor_indices)[0])
    return neighbor_indices
    

def generate_sphere_points(n):
    """
    Returns list of 3d coordinates of points on a sphere using the
    Golden Section Spiral algorithm.
    """
    points = []
    inc = math.pi * (3 - math.sqrt(5))
    offset = 2 / float(n)
    for k in range(int(n)):
        y = k * offset - 1 + (offset / 2)
        r = math.sqrt(1 - y*y)
        phi = k * inc
        points.append([math.cos(phi)*r, y, math.sin(phi)*r])
    return points
'''

#Original function
@profile
def find_neighbor_indices(atoms, probe, k):
    """
    Returns list of indices of atoms within probe distance to atom k. 
    """
    neighbor_indices = []
    atom_k = atoms[k]
    radius = atom_k.radius + probe + probe
    indices = range(k)
    indices.extend(range(k+1, len(atoms)))
    for i in indices:
        atom_i = atoms[i]
        dist = pos_distance(atom_k.pos, atom_i.pos)
        if dist < radius + atom_i.radius:
            neighbor_indices.append(i)
    return neighbor_indices
'''

#@profile
def calculate_asa(atoms, probe, n_sphere_point=960):
    """
    Returns list of accessible surface areas of the atoms, using the probe
    and atom radius to define the surface.
    """
    sphere_points = generate_sphere_points(n_sphere_point)
  
    points = np.array([ [a.pos.x, a.pos.y, a.pos.z] for a in atoms ]) #, dtype=np.float16) 
    radii = np.array([a.radius for a in atoms]) #, dtype=np.float16)

    const = 4.0 * math.pi / len(sphere_points)

    num_atoms = len(atoms)

    areas = np.zeros(num_atoms)

    for i in xrange(0, num_atoms):

        neighbor_indices = find_neighbor_indices(points, radii, probe, i)

        radius = probe + radii[i]

        point_i = points[i,:]

        n_accessible_point = 0

        for point in sphere_points:
            test_point = np.array(point)*radius + point_i
            neighbor_points = points[neighbor_indices, :]
            neighbor_radii_sq = (radii[neighbor_indices] + probe)**2

            diff_sq = np.sum((( neighbor_points - test_point)** 2), axis=1)
            dist_test = (diff_sq< neighbor_radii_sq)
            if dist_test.sum()==0:
                n_accessible_point += 1

        areas[i] = const*n_accessible_point*radius*radius 

    return areas


def main():
  import sys
  import getopt
  import molecule
  

  usage = \
  """  
  Calculates the total Accessible Surface Area (ASA) of atoms in a 
  PDB file. 

  Usage: asa.py -s n_sphere in_pdb [out_pdb]
  
  - out_pdb    PDB file in which the atomic ASA values are written 
               to the b-factor column.
               
  -s n_sphere  number of points used in generating the spherical
               dot-density for the calculation (default=960). The 
               more points, the more accurate (but slower) the 
               calculation. 

  """

  opts, args = getopt.getopt(sys.argv[1:], "n:")
  if len(args) < 1:
    print usage
    return
    
  mol = molecule.Molecule(args[0])
  atoms = mol.atoms()
  molecule.add_radii(atoms)

  n_sphere = 960 
  for o, a in opts:
    if '-n' in o:
      n_sphere = int(a)
      print "Points on sphere: ", n_sphere
  asas = calculate_asa(atoms, 1.4, n_sphere)
  print "%.1f angstrom squared." % sum(asas)

  if len(args) > 1:
    for asa, atom in zip(asas, atoms):
      atom.bfactor = asa
    mol.write_pdb(args[1])
  
  
if __name__ == "__main__":
  main()


    
