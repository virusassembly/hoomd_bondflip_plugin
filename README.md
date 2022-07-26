# Plugin for hoomd-blue
Folder 'bondflip_plugin' is a plugin for hoomd-blue. 
For a general plugin template, please refer to "example_plugin" in hoomd source: [**HOOMD-blue**
website](https://github.com/glotzerlab/hoomd-blue).


To use the plugin, put the folder 'bondflip_plugin' at the same path as 'example_plugin' and create a symbolic link to the plugin directory inside the hoomd subdirectory:

    $ cd /path/to/hoomd-blue/hoomd
    $ ln -s ../bondflip_plugin/bondflip_plugin bondflip_plugin
    $ cd ../build && make install

For details, please check [Plugins in HOOMD-Blue](https://hoomd-blue.readthedocs.io/en/v2.9.7/developer.html).

Next, compile hoomd following its instruction: 
[Compiling from source](INSTALLING.rst), notice that the version used here is [**2.9.7**](https://github.com/glotzerlab/hoomd-blue/releases/tag/v2.9.7).

Then the plugin is ready to be used.

Example for plugin bondflip:

    bondflip = hoomd.bondflip_plugin.update.bondflip(
        bond_neigh_bonds, 
        bond_dihedrals, 
        bond_faces,
        hoomd.group.type('M'), 
        logger, 
        period=100)

Note, bond_neigh_bonds, bond_dihedrals and bond_faces are generated in python script and pass to hoomd as arguments. 

Potential parameters are set through

    bondflip.set_params(bond_energy_name = 'harmonic',
        k = 20., 
        r0 = 1., 
        Da = 1.,
        Db = 1.,
        dihedral = True, 
        angle = False,
        kappa = 20., 
        phi_0 = np.pi, 
        T = 1.)

Two bond energy types: "harmonic" and "fene" is implemented in the plugin, please go to hoomd potential to check the parameter usage.
[bond potentials](https://hoomd-blue.readthedocs.io/en/v2.9.7/module-md-bond.html).


# Python script
An spherical vesicle example is included in the python_example folder, run python script by

    $ python main.py paramfile.json

Note, all parameters are included in "paramfile.json", and some alternative input files are provided in the folder.

