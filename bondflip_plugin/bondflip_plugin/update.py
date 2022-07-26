import hoomd
from hoomd import _hoomd
from hoomd.md import _md
from hoomd.bondflip_plugin import _bondflip_plugin
from hoomd.md import nlist as nl
import copy;



class bondflip(hoomd.update._updater):
    # Initialize the bondflip plugin
    #
    # \param bond_neigh_bonds : Dictionary of bond index and its four neighbor edges' indices.
    #        bond_dihedral: Dictionary of bond index and its dihedral index.
    #        bond_faces: Dictionary of bond index and its triangle indices.
    #        group : Group used to calculate the temperature.
    #        logger : Logger used to record the number of accepted bonds.
    #        period: Time steps between two bondflip attempts.
    def __init__(self, bond_neigh_bonds, bond_dihedral, bond_faces, group, logger, period=1):
        '''
        Creates an instance of the bondflip updater.
        :param bond_neigh_bonds: Dictionary of bond index and its four neighbor edges' indices.
        :param bond_dihedral: Dictionary of bond index and its dihedral index.
        :param bond_faces: Dictionary of bond index and its two triangle indices.
        :param group: Group used to calculate the temperature.
        :param logger: Logger used to record the number of accepted bonds.
        :param period: Time steps between two bondflip attempts.
        '''
        print("initializing bondflip in python")

        hoomd.util.print_status_line();

        # initialize base class
        hoomd.update._updater.__init__(self);

        self.bond_neigh_bonds = bond_neigh_bonds
        self.bond_dihedral = bond_dihedral
        self.bond_faces = bond_faces
        self.group = group
        self.logger = logger

        # create the compute thermo
        if group is hoomd.context.current.group_all:
            group_copy = copy.copy(group);
            group_copy.name = "__nvt_all";
            hoomd.util.quiet_status();
            self.thermo = hoomd.compute.thermo(group_copy);
            self.thermo.cpp_compute.setLoggingEnabled(False);
            hoomd.util.unquiet_status();
        else:
            self.thermo = hoomd.compute._get_unique_thermo(group=group);

        self.cpp_updater = _bondflip_plugin.BondFlipUpdater(hoomd.context.current.system_definition, 
            self.thermo.cpp_compute, 
            self.logger.cpp_analyzer, 
            self.bond_neigh_bonds, 
            self.bond_dihedral, 
            self.bond_faces);
        self.setupUpdater(period);

    def set_params(self,
                   bond_energy_name = "harmonic",
                   k = 1.0,
                   r0 = 1.0,
                   sigma = 0.0,
                   epsilon = 0.0,
                   Da = 1.0,
                   Db = 1.0,
                   dihedral = False,
                   angle = False,
                   kappa = 1.0,
                   phi_0 = 0.0,
                   T = 1.0):
        '''
        Sets parameters to the bondflip updater. 
        '''
        self.cpp_updater.set_params(bond_energy_name, k, r0, sigma, epsilon, Da, Db, dihedral, angle, kappa, phi_0, T)
    
    def get_bond_neigh_bonds(self):
        return self.cpp_updater.get_bond_neigh_bonds()
    
    def get_bond_dihedrals(self):
        return self.cpp_updater.get_bond_dihedrals()
    
    def get_bond_faces(self):
        return self.cpp_updater.get_bond_faces()

