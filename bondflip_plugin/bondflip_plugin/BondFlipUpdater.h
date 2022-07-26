#ifndef _BONDFLIP_UPDATER_H_
#define _BONDFLIP_UPDATER_H_

/*! \file BondFlipUpdater.h
*/
#include "hoomd/ComputeThermo.h"
#include "hoomd/Logger.h"
#include "hoomd/ParticleGroup.h"
#include <hoomd/Updater.h>
#include <map>
#include <random>
#ifndef NVCC
#include <hoomd/extern/pybind/include/pybind11/pybind11.h>
#include <hoomd/extern/pybind/include/pybind11/stl.h>

#endif

class BondFlipUpdater : public Updater
    {
    public:
        //! Constructor
        BondFlipUpdater(std::shared_ptr<SystemDefinition> sysdef,
                std::shared_ptr<ComputeThermo> thermo,
                std::shared_ptr<Logger> logger,
                std::map<unsigned int, std::vector<unsigned int>> bond_neigh_bonds,
                std::map<unsigned int, unsigned int> bond_dihedral,
                std::map<unsigned int, std::vector<unsigned int>> bond_faces);
        void set_params(std::string bond_energy_name,
                        Scalar k,
                        Scalar r0,
                        Scalar sigma,
                        Scalar epsilon,
                        Scalar diameter_a,
                        Scalar diameter_b,
                        bool dihedral,
                        bool angle,
                        Scalar kappa,
                        Scalar phi_0,
                        Scalar T);
        //! Take one timestep forward
        virtual void update(unsigned int timestep);
        //! Returns a list of log quantities this compute calculates
        virtual std::vector< std::string > getProvidedLogQuantities(void);
        //! Calculates the requested log value and returns it
        virtual Scalar getLogValue(const std::string& quantity, unsigned int timestep);
        std::map<unsigned int, std::vector<unsigned int> > get_bond_neigh_bonds();
        std::map<unsigned int, unsigned int> get_bond_dihedrals();
        std::map<unsigned int, std::vector<unsigned int> > get_bond_faces();
        
        std::map<unsigned int, pybind11::object> m_callback;

    protected:
        std::mt19937 m_rnd;
        unsigned int m_seed;                               
        std::string m_bond_energy_name;
        Scalar m_k = 1.0;
        Scalar m_r0 = 1.0;
        Scalar m_sigma = 0.0;
        Scalar m_epsilon = 0.0;
        Scalar m_diameter_a = 1.0;
        Scalar m_diameter_b = 1.0;
        Scalar m_T = 1.0;
        bool m_dihedral = false;
        bool m_angle = false;
        Scalar m_kappa = 1.0;
        Scalar m_phi_0 = 3.1415926;
        int m_accepted = 0;

        bool is_connected(unsigned int, unsigned int);
        void find_opposite_vertex(unsigned int, unsigned int, unsigned int, unsigned int, unsigned int &, unsigned int &);
        void flip_bond(unsigned int, Scalar);
        void printcheck(unsigned int);
        
        Scalar distance(BoxDim, ArrayHandle<Scalar4>, unsigned int, unsigned int);
        bool face_are_neighbor(unsigned int f0, unsigned int f1);
        void update_neigh_bonds(unsigned int, unsigned int, std::vector<unsigned int>, std::map<unsigned int, unsigned int>);
        bool get_bond_energy(BoxDim box, ArrayHandle<Scalar4> h_pos, unsigned int idx_a, unsigned int idx_b, Scalar& bond_eng);
        Scalar get_dihedral_energy(BoxDim, ArrayHandle<Scalar4>, unsigned int idx_a, unsigned int idx_b, unsigned int idx_c, unsigned int idx_d);

        std::map<unsigned int, std::vector<unsigned int> > m_bond_neigh_bonds; //!< Dictionary of bond index and four neighbor edges' indices
        std::map<unsigned int, unsigned int> m_bond_dihedral; //!< Dictionary of bond index and its dihedral index, 32767 if bond is on boundary
        std::map<unsigned int, std::vector<unsigned int>> m_bond_faces;

        std::map<unsigned int, unsigned int> m_vertex_bond_num; //!< Dictionary of vertex index and number of connected bonds
        std::shared_ptr<BondData> m_bonds;
        std::shared_ptr<AngleData> m_angles;
        std::shared_ptr<DihedralData> m_dihedrals;
        const std::shared_ptr<ComputeThermo> m_thermo;
        const std::shared_ptr<Logger> m_logger; 
        std::vector<std::string> m_loggable_quantities;
    };

//! Export the BondFlipUpdater class to python
void export_BondFlipUpdater(pybind11::module& m);

#ifdef ENABLE_CUDA

//! A GPU accelerated version of the BondFlipUpdater.
class BondFlipUpdaterGPU : public BondFlipUpdater
    {
    public:
        //! Constructor
        BondFlipUpdaterGPU(std::shared_ptr<SystemDefinition> sysdef,
                    std::shared_ptr<ComputeThermo> thermo,
                    std::shared_ptr<Logger> logger,
                    std::map<unsigned int, std::vector<unsigned int>> bond_neigh_bonds,
                    std::map<unsigned int, unsigned int> bond_dihedral,
                    std::map<unsigned int, std::vector<unsigned int>> bond_faces);
        //! Take one timestep forward
        virtual void update(unsigned int timestep);
    };

//! Export the BondFlipUpdaterGPU class to python
void export_BondFlipUpdaterGPU(pybind11::module& m);

#endif // ENABLE_CUDA

#endif // _BONDFLIP_UPDATER_H_