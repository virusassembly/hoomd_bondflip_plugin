# @Author: Siyu Li(sli032@ucr.edu)
# @Date:   2022-06-07 15:39:41
# @Last Modified by:   Siyu Li
# @Last Modified time: 2022-06-07 15:42:45
import numpy as np

def assign_snap_particles(snap, obj, typedic, cover):
	# particles assign in range:[pstart,pend)
	pstart = snap.particles.N
	pend = pstart + obj.n
	snap.particles.resize(pend)
	print("pstart, pend:",pstart,pend)

	for i, (cor, d, typename) in enumerate(zip(obj.particles, obj.dlist, obj.ptypelist), start = pstart):
		# print('cor:{},pid:{}'.format(cor,pid))
		snap.particles.diameter[i] = d
		snap.particles.position[i] = cor
		snap.particles.typeid[i] = typedic[typename]
		snap.particles.body[i] = obj.bodyid
		
	return pstart, pend

def assign_snap_bonds(snap, obj, typedic, pstart, cover):
	bstart=snap.bonds.N
	snap.bonds.resize(bstart + len(obj.bonds))
	print('snap.bonds:', snap.bonds.types)
	
	for i, bond in enumerate(obj.bonds, start = bstart):
		snap.bonds.group[i] = pstart + bond
		snap.bonds.typeid[i] = typedic[obj.btypename]

def assign_snap_angles(snap, obj, typedic, pstart, cover):
	dstart = snap.angles.N
	snap.angles.resize(dstart + len(obj.angles))
	
	for i, angle in enumerate(obj.angles, start = dstart):
		snap.angles.group[i] = pstart + angle
		snap.angles.typeid[i] = typedic[obj.atypename]

def assign_snap_dihedrals(snap, obj, typedic, pstart, cover):
	dstart = snap.dihedrals.N
	snap.dihedrals.resize(dstart + len(obj.dihedrals))
	
	for i, dihedral in enumerate(obj.dihedrals, start = dstart):
		snap.dihedrals.group[i] = pstart + dihedral
		snap.dihedrals.typeid[i] = typedic[obj.dtypename]

def update_snap(snap, obj, typedic, bond = True, angle = False, dihedral = True):
	pstart,pend=assign_snap_particles(snap, obj,typedic, False)
	if bond:
		assign_snap_bonds(snap,obj,typedic,pstart,False)
	if angle:
		assign_snap_angles(snap,obj,typedic,pstart,False)
	if dihedral:
		assign_snap_dihedrals(snap,obj,typedic,pstart,False)
	return pstart, pend