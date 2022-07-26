# @Author: Siyu Li(sli032@ucr.edu)
# @Date:   2022-06-07 15:39:41
# @Last Modified by:   Siyu Li
# @Last Modified time: 2022-06-07 15:42:52
import sys
import hoomd
import hoomd.md
import hoomd.bondflip_plugin
import geometry as gm
import system
import numpy as np
import os
import hoomd.deprecated as hdepr
import json, glob


def snapprint(snap):
    np.set_printoptions(threshold=np.inf)
    print('particle position:\n{} \
            \nparticle typeid:\n{} \
            \npartical body:\n{} \
            \nparticale diameter:\n{}' \
            .format(np.array(snap.particles.position[:]),np.array(snap.particles.typeid[:]),np.array(snap.particles.body[:]),np.array(snap.particles.diameter[:])))
    print('bond groups:\n{} \
            \nbond typeid:\n{}' \
            .format(np.array(snap.bonds.group[:]),np.array(snap.bonds.typeid[:])))
    print('angle groups:\n{} \
            \nangle typeid:\n{}' \
            .format(np.array(snap.angles.group[:]),np.array(snap.angles.typeid[:])))
    print('dihedral groups:\n{} \
            \ndihedral typeid:\n{}' \
            .format(np.array(snap.dihedrals.group[:]),np.array(snap.dihedrals.typeid[:])))

def set_system(mem):   
    logquantities = ['temperature', 'potential_energy',
    'translational_kinetic_energy', 'rotational_kinetic_energy',
    'bond_harmonic_energy', 'dihedral_harmonic_energy',
    'accepted_bonds']
    logger = hoomd.analyze.log(filename="energy.log",
        quantities=logquantities,
        period=param['record_period'], overwrite=False)

    bond_harmonic = hoomd.md.bond.harmonic()
    dihedral_force = hoomd.md.dihedral.harmonic()
    bond_harmonic.bond_coeff.set('mbond', k = param['mbondk'], r0 = param['mbond_r0'])    
    dihedral_force.dihedral_coeff.set('mdihedral', k = 2 * param['mdihedk'], d = -1, n = 1, phi_0 = np.pi - param['mtheta0'])

    bondflip = hoomd.bondflip_plugin.update.bondflip(mem.bond_neigh_bonds, 
        mem.bond_dihedrals, 
        mem.bond_faces,
        hoomd.group.all(), 
        logger, 
        period=param['bondflip_period'])
    
    bondflip.set_params(bond_energy_name = 'harmonic',
        k = param['mbondk'], 
        r0 = param['mbond_r0'], 
        Da = 1.,
        Db = 1.,
        dihedral = True, 
        angle = False,
        kappa = param['mdihedk'], 
        phi_0 = np.pi - param['mtheta0'], 
        T = param['bondflip_kT'])
        

def main():
    hoomd.context.initialize("--mode=cpu")

    # generate simulation box
    snap=hoomd.data.make_snapshot(
        N = 0,
        box = hoomd.data.boxdim(
            Lx = param['boxL'],
            Ly = param['boxL'],
            Lz = param['boxL']))
    snap.particles.types = ['M']
    snap.bonds.types = ['mbond']
    snap.angles.types = ['mangle']
    snap.dihedrals.types = ['mdihedral']
    typedic = {'M':0, 'mbond':0, 'mangle':0, 'mdihedral':0}

    # generate membrane
    mem = gm.membrane.mem()
    mem.frame_in(param['memfile'])
    mem.update_mem_info()
    mem.update_bond_flip_info()
    system.snap.update_snap(snap, mem, typedic, bond = True, angle = False, dihedral = True)

    # snapprint(snap)

    # generate hoomd system
    hoomd.init.read_snapshot(snap)
    
    set_system(mem)

    hoomd.md.integrate.mode_standard(dt=param['dt'])
    integrator = hoomd.md.integrate.langevin(group = hoomd.group.type('M'), kT = param['kT'], seed = param['rseed'])

    hoomd.dump.gsd(filename="md.gsd", group=hoomd.group.all(), period=None, static=[])
    for _ in range(param['runtime']):
        hoomd.run(param['record_period'])
        hoomd.dump.gsd(filename="md.gsd", group=hoomd.group.all(), period=None, static=[])

param=sys.argv[1]
print(param)
with open(param) as f: param=json.load(f)
main()