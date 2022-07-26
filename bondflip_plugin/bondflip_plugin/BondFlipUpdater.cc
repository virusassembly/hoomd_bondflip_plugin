#include "BondFlipUpdater.h"
#ifdef ENABLE_CUDA
#include "BondFlipUpdater.cuh"
#endif
#include "hoomd/Saru.h"
#include <math.h>

#define ONEDGE 32767

// windows defines a macro min and max
// #undef min
// #undef max

using namespace hoomd;
/*! \file BondFlipUpdater.cc
    \BondFlip Updater clip the shared edge of two adjacent triangles and reconnect to the diagonal vertices with a probability defined by the boltzmann factor.
*/


/*! \param sysdef: System to perform bond flip
    \param thermo: ComputeThermo object to get system temperature
*/
BondFlipUpdater::BondFlipUpdater(std::shared_ptr<SystemDefinition> sysdef, 
            std::shared_ptr<ComputeThermo> thermo,
            std::shared_ptr<Logger> logger,std::map<unsigned int, std::vector<unsigned int>> bond_neigh_bonds,
            std::map<unsigned int, unsigned int> bond_dihedral,
            std::map<unsigned int, std::vector<unsigned int>> bond_faces)
        : Updater(sysdef), m_thermo(thermo), m_logger(logger)
    {
        printf("Initializing BondFlipUpdater.....");
        m_bond_neigh_bonds = bond_neigh_bonds; // ..
        m_bond_dihedral = bond_dihedral;
        m_bond_faces = bond_faces;
        // get system bonds data
        m_bonds = m_sysdef->getBondData();
        m_angles = m_sysdef->getAngleData();
        m_dihedrals = m_sysdef->getDihedralData();
        // random number generator
        m_seed = 1234;
        m_rnd.seed(std::mt19937::result_type(m_seed));

        unsigned int bondN = m_bonds->getN();
        // loop all bonds to get m_vertex_bond_num
        for (unsigned int i = 0; i < bondN; i++)
        {
            Bond bond = m_bonds -> getGroupByTag(i);
            unsigned int a = bond.a;
            unsigned int b = bond.b;
            m_vertex_bond_num[a]++;
            m_vertex_bond_num[b]++;
        }
    }


void BondFlipUpdater::set_params(std::string bond_energy_name,
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
                               Scalar T)
    {
        m_bond_energy_name = bond_energy_name;
        m_k = k;
        m_r0 = r0;
        m_sigma = sigma;
        m_epsilon = epsilon;
        m_diameter_a = diameter_a;
        m_diameter_b = diameter_b;
        m_dihedral = dihedral;
        m_angle = angle;
        m_kappa = kappa;
        m_phi_0 = phi_0;
        m_T = T;
    }

bool BondFlipUpdater::is_connected(unsigned int b0, unsigned int b1)
{
    Bond bond0 = m_bonds -> getGroupByTag(b0);
    Bond bond1 = m_bonds -> getGroupByTag(b1);
    return bond0.a == bond1.a || bond0.a == bond1.b || bond0.b == bond1.a || bond0.b == bond1.b;
}
bool BondFlipUpdater::face_are_neighbor(unsigned int f0, unsigned int f1)
{
    Angle angle0 = m_angles -> getGroupByTag(f0);
    Angle angle1 = m_angles -> getGroupByTag(f1);
    std::set<unsigned int> v = {angle0.a, angle0.b, angle0.c, angle1.a, angle1.b, angle1.c};
    return v.size() == 4;
}
void BondFlipUpdater::find_opposite_vertex(unsigned int bi, unsigned int bi_op, unsigned int a, unsigned int b, unsigned int &c, unsigned int &d)
{
    Bond bond = m_bonds -> getGroupByTag(bi);
    Bond bond_op = m_bonds -> getGroupByTag(bi_op);
    unsigned int tot_v = bond.a + bond.b + bond_op.a + bond_op.b;
    if (bond.a == a || bond.b == a)
        c = bond.a + bond.b - a;
    else
        c = bond.a + bond.b - b;
    d = tot_v - a - b - c;
}
Scalar BondFlipUpdater::distance(BoxDim box, ArrayHandle<Scalar4> h_pos, unsigned int idx_a, unsigned int idx_b)
{
    vec3<Scalar> pi(h_pos.data[idx_a].x, h_pos.data[idx_a].y, h_pos.data[idx_a].z);
    vec3<Scalar> pj(h_pos.data[idx_b].x, h_pos.data[idx_b].y, h_pos.data[idx_b].z);
    vec3<Scalar> drScalar(pj - pi);
    drScalar = vec3<Scalar>(box.minImage(vec_to_scalar3(drScalar)));
    Scalar dist = sqrt(drScalar.x*drScalar.x + drScalar.y*drScalar.y + drScalar.z*drScalar.z);
    return dist;
}
void BondFlipUpdater::update_neigh_bonds(unsigned int curr_b, unsigned int new_b, std::vector<unsigned int> neigh_bonds, std::map<unsigned int, unsigned int> op_bond_dic)
{
    if (new_b != curr_b)
    {
        // update m_bond_neigh_bonds for curr_b
        auto entry = m_bond_neigh_bonds.find(curr_b);
        if (entry != end(m_bond_neigh_bonds))
        {
            auto const value = std::move(entry->second);
            m_bond_neigh_bonds.erase(entry);
            m_bond_neigh_bonds.insert({new_b, std::move(value)});
        }
    }
    
    // update m_bond_neigh_bonds for neighbor bonds
    for (unsigned int bi: neigh_bonds)
        if (m_bond_neigh_bonds.find(bi) != m_bond_neigh_bonds.end())
        {
            for(unsigned int k = 0; k < 4; k++)
            {
                if (m_bond_neigh_bonds[bi][k] == curr_b)
                    m_bond_neigh_bonds[bi][k] = new_b;
                
                if (op_bond_dic.find(m_bond_neigh_bonds[bi][k]) != op_bond_dic.end())
                {
                    m_bond_neigh_bonds[bi][k] = op_bond_dic[m_bond_neigh_bonds[bi][k]];
                }
            }
        }
}

Scalar BondFlipUpdater::get_dihedral_energy(BoxDim box, ArrayHandle<Scalar4> h_pos, unsigned int idx_a, unsigned int idx_b, unsigned int idx_c, unsigned int idx_d)
{
    // printf("calculating dihedral energy...");
    // Adopt from HarmonicDihedralForceCompute.cc
    Scalar3 dab;
    dab.x = h_pos.data[idx_a].x - h_pos.data[idx_b].x;
    dab.y = h_pos.data[idx_a].y - h_pos.data[idx_b].y;
    dab.z = h_pos.data[idx_a].z - h_pos.data[idx_b].z;

    Scalar3 dcb;
    dcb.x = h_pos.data[idx_c].x - h_pos.data[idx_b].x;
    dcb.y = h_pos.data[idx_c].y - h_pos.data[idx_b].y;
    dcb.z = h_pos.data[idx_c].z - h_pos.data[idx_b].z;

    Scalar3 ddc;
    ddc.x = h_pos.data[idx_d].x - h_pos.data[idx_c].x;
    ddc.y = h_pos.data[idx_d].y - h_pos.data[idx_c].y;
    ddc.z = h_pos.data[idx_d].z - h_pos.data[idx_c].z;

    // apply periodic boundary conditions
    dab = box.minImage(dab);
    dcb = box.minImage(dcb);
    ddc = box.minImage(ddc);

    Scalar3 dcbm;
    dcbm.x = -dcb.x;
    dcbm.y = -dcb.y;
    dcbm.z = -dcb.z;

    dcbm = box.minImage(dcbm);

    Scalar aax = dab.y*dcbm.z - dab.z*dcbm.y;
    Scalar aay = dab.z*dcbm.x - dab.x*dcbm.z;
    Scalar aaz = dab.x*dcbm.y - dab.y*dcbm.x;

    Scalar bbx = ddc.y*dcbm.z - ddc.z*dcbm.y;
    Scalar bby = ddc.z*dcbm.x - ddc.x*dcbm.z;
    Scalar bbz = ddc.x*dcbm.y - ddc.y*dcbm.x;

    Scalar raasq = aax*aax + aay*aay + aaz*aaz;
    Scalar rbbsq = bbx*bbx + bby*bby + bbz*bbz;
    Scalar rgsq = dcbm.x*dcbm.x + dcbm.y*dcbm.y + dcbm.z*dcbm.z;
    Scalar rg = sqrt(rgsq);

    Scalar rginv, raa2inv, rbb2inv;
    rginv = raa2inv = rbb2inv = Scalar(0.0);
    if (rg > Scalar(0.0)) rginv = Scalar(1.0)/rg;
    if (raasq > Scalar(0.0)) raa2inv = Scalar(1.0)/raasq;
    if (rbbsq > Scalar(0.0)) rbb2inv = Scalar(1.0)/rbbsq;
    Scalar rabinv = sqrt(raa2inv*rbb2inv);

    Scalar c_abcd = (aax*bbx + aay*bby + aaz*bbz)*rabinv;
    Scalar s_abcd = rg*rabinv*(aax*ddc.x + aay*ddc.y + aaz*ddc.z);

    if (c_abcd > 1.0) c_abcd = 1.0;
    if (c_abcd < -1.0) c_abcd = -1.0;

    // printf("cos: %f, sin:%f\n", c_abcd, s_abcd);

    Scalar sin_phi_0 = fast::sin(m_phi_0);
    Scalar cos_phi_0 = fast::cos(m_phi_0);
    Scalar p = Scalar(1.0) - c_abcd*cos_phi_0 - s_abcd*sin_phi_0;

    Scalar energy =  m_kappa * p;
    // printf("dihedral energy is %f\n", energy);
    return energy;
}

bool BondFlipUpdater::get_bond_energy(BoxDim box, ArrayHandle<Scalar4> h_pos, unsigned int idx_a, unsigned int idx_b, Scalar& bond_eng)
{
    // Adopt from EvaluatorBondFENE.h
    Scalar3 rab;
    rab.x = h_pos.data[idx_a].x - h_pos.data[idx_b].x;
    rab.y = h_pos.data[idx_a].y - h_pos.data[idx_b].y;
    rab.z = h_pos.data[idx_a].z - h_pos.data[idx_b].z;

    // apply periodic boundary conditions
    rab = box.minImage(rab);
    Scalar rsq = rab.x*rab.x + rab.y*rab.y + rab.z*rab.z;
    if (m_bond_energy_name.compare("harmonic") == 0)
    {
        Scalar r = sqrt(rsq);
        bond_eng = Scalar(0.5) * m_k * (r - m_r0) * (r - m_r0);
    }
    else if (m_bond_energy_name.compare("fene") == 0)
    {
        Scalar rtemp = sqrt(rsq) - m_diameter_a/2 - m_diameter_b/2 + Scalar(1.0);
        rsq = rtemp*rtemp;

        Scalar r2inv = Scalar(1.0)/rsq;
        Scalar r6inv = r2inv * r2inv * r2inv;
        Scalar sigma2 = m_sigma * m_sigma;
        Scalar sigma6 = sigma2 * sigma2 * sigma2;

        Scalar pair_eng = Scalar(0.0);

        if (sigma6 * r6inv > Scalar(0.5))
            {
            pair_eng = (4 * m_epsilon * sigma6 * r6inv * (sigma6*r6inv - 1) + m_epsilon);
            }
        
        if (rsq >= m_r0*m_r0) 
            return false;
        else
            bond_eng = -Scalar(0.5) * m_k * (m_r0 * m_r0) * log(Scalar(1.0) - rsq/(m_r0 * m_r0)) + pair_eng;
    }
    else
        bond_eng = 0.;
        // m_exec_conf->msg->warning() << "bond energy name not found!";
    return true;
}

//! Helper function to generate a [0..1] float
/*! \param rnd Random number generator to use
*/
static Scalar random01(std::mt19937& rnd)
    {
    unsigned int val = rnd();

    double ranval = ((double)val - (double)rnd.min()) / ( (double)rnd.max() - (double)rnd.min() );
    return Scalar(ranval);
    }
/*! Perform the bond flip
    \param i index of the chosen bond
*/
// Note, everything is dealt with the bond tag!
void BondFlipUpdater::flip_bond(unsigned int curr_b, Scalar curr_T)
{    
    // check if the bond is on boundary
    if (m_bond_neigh_bonds.find(curr_b) == m_bond_neigh_bonds.end())
    {
        // m_exec_conf->msg->warning() << "Bond is on Boundary!";
        return;
    }

    Bond bond = m_bonds -> getGroupByTag(curr_b);
    unsigned int a = bond.a;
    unsigned int b = bond.b;
    unsigned int c, d;
    unsigned int b_typeid = bond.type;
    std::vector<unsigned int> neigh_bonds = m_bond_neigh_bonds[curr_b];
    unsigned int b0 = neigh_bonds[0];
    unsigned int b1 = neigh_bonds[1];
    unsigned int b2 = neigh_bonds[2];
    unsigned int b3 = neigh_bonds[3];
    std::map<unsigned int, unsigned int> op_bond_dic;
    if (!is_connected(b0, b1))
    {
        op_bond_dic[b0] = b1;
        op_bond_dic[b1] = b0;
        op_bond_dic[b2] = b3;
        op_bond_dic[b3] = b2;
    }
    else if (!is_connected(b0, b2))
    {
        op_bond_dic[b0] = b2;
        op_bond_dic[b2] = b0;
        op_bond_dic[b1] = b3;
        op_bond_dic[b3] = b1;
    }
    else if (!is_connected(b0, b3))
    {
        op_bond_dic[b0] = b3;
        op_bond_dic[b3] = b0;
        op_bond_dic[b1] = b2;
        op_bond_dic[b2] = b1;
    }
    else
    {
        m_exec_conf->msg->warning() << "Cannot find opposite bonds!";
    }

    find_opposite_vertex(b0, op_bond_dic[b0], a, b, c, d);
    // printf("original vertex: %d, %d, opposite vertex: %d, %d\n",a, b, c, d);

    ArrayHandle<Scalar4> h_pos(m_pdata->getPositions(), access_location::host, access_mode::read);
    ArrayHandle< unsigned int > h_rtag(m_pdata->getRTags(), access_location::host, access_mode::read);
    BoxDim box = m_pdata->getBox();
    Scalar E_old = Scalar(0.0);
    Scalar E_new = Scalar(0.0);
    bool eva_old = get_bond_energy(box, h_pos, h_rtag.data[a], h_rtag.data[b], E_old);
    bool eva_new = get_bond_energy(box, h_pos, h_rtag.data[c], h_rtag.data[d], E_new);
    

    if (eva_old && eva_new)
    {
        // printf("eva_old:%d, eva_new:%d", eva_old, eva_new);
        
        if (m_dihedral)
        {
            unsigned int curr_d = m_bond_dihedral[curr_b];
            Dihedral diheral = m_dihedrals -> getGroupByTag(curr_d);
            unsigned int v0 = diheral.a;
            unsigned int v1 = diheral.b;
            unsigned int v2 = diheral.c;
            unsigned int v3 = diheral.d;
            unsigned int d_typeid = diheral.type;
            E_old += get_dihedral_energy(box, h_pos, h_rtag.data[v0], h_rtag.data[v1], h_rtag.data[v2], h_rtag.data[v3]);
            E_new += get_dihedral_energy(box, h_pos, h_rtag.data[v1], h_rtag.data[v3], h_rtag.data[v0], h_rtag.data[v2]);
            // printf("E_old: %f, E_new: %f",E_old,E_new);
            // printf("v0:%d, v1:%d, v2:%d, v3:%d, v0 cor:%f", v0, v1, v2, v3, h_pos.data[h_rtag.data[v0]].x);
            for (auto bi: neigh_bonds)
            {
                if (m_bond_dihedral.find(bi) != m_bond_dihedral.end())
                {
                    unsigned int di = m_bond_dihedral[bi];
                    Dihedral diheral = m_dihedrals -> getGroupByTag(di);
                    unsigned int d0 = diheral.a;
                    unsigned int d1 = diheral.b;
                    unsigned int d2 = diheral.c;
                    unsigned int d3 = diheral.d;
                    E_old += get_dihedral_energy(box, h_pos, h_rtag.data[d0], h_rtag.data[d1], h_rtag.data[d2], h_rtag.data[d3]);
                    if (d0 == a || d0 == b)
                    {
                        unsigned int other_v = v0 + v1 + v2 + v3 - d0 - d1 - d2;
                        E_new += get_dihedral_energy(box, h_pos, h_rtag.data[other_v], h_rtag.data[d1], h_rtag.data[d2], h_rtag.data[d3]);
                    }
                    else if (d3 == a || d3 == b)
                    {
                        unsigned int other_v = v0 + v1 + v2 + v3 - d1 - d2 - d3;
                        E_new += get_dihedral_energy(box, h_pos, h_rtag.data[d0], h_rtag.data[d1], h_rtag.data[d2], h_rtag.data[other_v]);
                    }
                    else
                    {
                        m_exec_conf->msg->warning() << "dihedral order seems wrong!";
                    }
                }
            }
        }
        
        Scalar dE = E_new - E_old;
        Scalar mb_stats = exp(-dE/m_T);
        Scalar temp_rand = random01(m_rnd);

        // printf("E_old:%f, E_new:%f, dE:%f, prob:%f", E_old, E_new, dE, mb_stats);

        if (m_vertex_bond_num[a] > 4 && m_vertex_bond_num[b] > 4 && m_vertex_bond_num[c] < 7 && m_vertex_bond_num[d] < 7)
            if (mb_stats > temp_rand)
            {
                // printf("bondflip accepted");
                // std::cout<<"before"<<std::endl;
                // printcheck(curr_b);
                m_bonds -> removeBondedGroup(curr_b);
                unsigned int new_b = m_bonds -> addBondedGroup(Bond(b_typeid, c, d));
                // printf("old bond %d is attempted to flip to new bond %d.\n", curr_b, new_b);
                m_vertex_bond_num[a]--;
                m_vertex_bond_num[b]--;
                m_vertex_bond_num[c]++;
                m_vertex_bond_num[d]++;

                update_neigh_bonds(curr_b, new_b, neigh_bonds, op_bond_dic);

                if (m_dihedral)
                {
                    unsigned int curr_d = m_bond_dihedral[curr_b];
                    Dihedral diheral = m_dihedrals -> getGroupByTag(curr_d);
                    unsigned int v0 = diheral.a;
                    unsigned int v1 = diheral.b;
                    unsigned int v2 = diheral.c;
                    unsigned int v3 = diheral.d;
                    unsigned int d_typeid = diheral.type;
                    m_dihedrals -> removeBondedGroup(curr_d);
                    unsigned int new_d = m_dihedrals -> addBondedGroup(Dihedral(d_typeid, v1, v3, v0, v2));
                    assert(curr_d == new_d);
                    for (auto bi: neigh_bonds)
                    {
                        if (m_bond_dihedral.find(bi) != m_bond_dihedral.end())
                        {
                            unsigned int di = m_bond_dihedral[bi];
                            Dihedral diheral = m_dihedrals -> getGroupByTag(di);
                            unsigned int d0 = diheral.a;
                            unsigned int d1 = diheral.b;
                            unsigned int d2 = diheral.c;
                            unsigned int d3 = diheral.d;
                            if (d0 == a || d0 == b)
                            {
                                unsigned int new_v = v0 + v1 + v2 + v3 - d0 - d1 - d2;
                                m_dihedrals -> removeBondedGroup(di);
                                new_d = m_dihedrals -> addBondedGroup(Dihedral(d_typeid, new_v, d1, d2, d3));
                            }
                            else if (d3 == a || d3 == b)
                            {
                                unsigned int new_v = v0 + v1 + v2 + v3 - d1 - d2 - d3;
                                m_dihedrals -> removeBondedGroup(di);
                                new_d = m_dihedrals -> addBondedGroup(Dihedral(d_typeid, d0, d1, d2, new_v));
                            }
                            assert(di == new_d);
                        }
                    }
                }

                if (m_angle)
                {
                    unsigned int face0 = m_bond_faces[curr_b][0];
                    unsigned int face1 = m_bond_faces[curr_b][1];
                    Angle angle0 = m_angles -> getGroupByTag(face0);
                    Angle angle1 = m_angles -> getGroupByTag(face1);

                    unsigned int an_typeid0 = angle0.type;
                    unsigned int an_typeid1 = angle1.type;
                    m_angles -> removeBondedGroup(face0);
                    unsigned int new_face0 = m_angles -> addBondedGroup(Angle(an_typeid0, a, c, d));
                    assert(face0 == new_face0);
                    m_angles -> removeBondedGroup(face1);
                    unsigned int new_face1 = m_angles -> addBondedGroup(Angle(an_typeid1, b, c, d));
                    assert(face1 == new_face1);

                    for (auto bi: neigh_bonds)
                    {
                        if (m_bond_faces.find(bi) != m_bond_faces.end())
                        {
                            unsigned int f0 = m_bond_faces[bi][0];
                            unsigned int f1 = m_bond_faces[bi][1];
                            if (f0 == face0 || f0 == face1)
                            {
                                if (face_are_neighbor(f1, face0))
                                    m_bond_faces[bi][0] = face0;
                                else
                                    m_bond_faces[bi][0] = face1;
                            }
                            else if (f1 == face0 || f1 == face1)
                            {
                                if (face_are_neighbor(f0, face0))
                                    m_bond_faces[bi][1] = face0;
                                else
                                    m_bond_faces[bi][1] = face1;
                            }      
                            else
                                m_exec_conf->msg->warning() << "bond_faces of neigh bonds are wrong!";   
                        }             
                    }
                }
                // std::cout<<"after"<<std::endl;
                // printcheck(curr_b);


                m_accepted++;
                // printf("bond flip accepted!");

                // printf("manual cal E_old: %f, E_new: %f, dE: %f, mb_stats: %f rand %f m_T:%f\n",E_old,E_new,dE,mb_stats,temp_rand,m_T);
            }
        else
        {
            // printf("bondflip refuse");
        }
    }
            
}


/*! Perform N times bondflips, where N is the number of total bonds.
    \param timestep Current time step of the simulation
*/
void BondFlipUpdater::update(unsigned int timestep)
    {
        // m_accepted = 0;
        if (m_prof) m_prof->push("BondFlipUpdater");
        // compute the current thermodynamic properties and get the temperature
        m_thermo->compute(timestep);
        Scalar curr_T = m_thermo->getTranslationalTemperature();
        
        unsigned int bondN = m_bonds->getN();
        for (unsigned int i = 0; i < bondN; i++)
        {
           // random pick on a bond
           unsigned int btag = round(random01(m_rnd)*(bondN - 1));
           // printf("try to flip bond %d.\n",btag);
           flip_bond(btag, curr_T);
        }
        if (m_prof) m_prof->pop();
    }

void BondFlipUpdater::printcheck(unsigned int btag)
{
    unsigned int bondN = m_bonds->getN();
    unsigned int angleN = m_angles->getN();
    unsigned int dihedralN = m_dihedrals->getN();
    ArrayHandle<unsigned int> bond_tags(m_bonds->getTags(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> bond_rtags(m_bonds->getRTags(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> dihedral_rtags(m_dihedrals->getRTags(), access_location::host, access_mode::read);
    ArrayHandle<unsigned int> angle_rtags(m_angles->getRTags(), access_location::host, access_mode::read);
    // for (unsigned int bondi = 0; bondi< bondN; bondi++)
    {
        std::cout<<"bond tag "<<btag<< "reverse index " << bond_rtags.data[btag] <<std::endl;
        if (m_bond_neigh_bonds.find(btag) != m_bond_neigh_bonds.end())
        {
            std::cout<< "neigh bonds:"<< m_bond_neigh_bonds[btag][0]<< " " << m_bond_neigh_bonds[btag][1]<< " "<<m_bond_neigh_bonds[btag][2]<<" "<<m_bond_neigh_bonds[btag][3]<<std::endl;
            unsigned int dtag = m_bond_dihedral[btag];
            Dihedral diheral = m_dihedrals -> getGroupByTag(dtag);
            unsigned int v0 = diheral.a;
            unsigned int v1 = diheral.b;
            unsigned int v2 = diheral.c;
            unsigned int v3 = diheral.d;
            std::cout<<"dihedral tag "<<dtag<< "reverse index " << dihedral_rtags.data[dtag] <<" vertex:"<<v0<<","<<v1<<","<<v2<<","<<v3<<std::endl;
            unsigned int face0 = m_bond_faces[btag][0];
            unsigned int face1 = m_bond_faces[btag][1];
            Angle angle0 = m_angles -> getGroupByTag(face0);
            Angle angle1 = m_angles -> getGroupByTag(face1);
            std::cout<<"face0 tag "<<face0<< "reverse index " << angle_rtags.data[face0] <<" vertex:"<<angle0.a<<","<<angle0.b<<","<<angle0.c<<std::endl;
            std::cout<<"face1 tag "<<face1<< "reverse index " << angle_rtags.data[face1] <<" vertex:"<<angle1.a<<","<<angle1.b<<","<<angle1.c<<std::endl;
        }
        
    }
    // for (unsigned int dihedrali = 0; dihedrali< dihedralN; dihedrali++)
    // {
    //     unsigned int dtag = dihedral_tags.data[dihedrali];
        
    //     Dihedral diheral = m_dihedrals -> getGroupByTag(dtag);
    //     unsigned int v0 = diheral.a;
    //     unsigned int v1 = diheral.b;
    //     unsigned int v2 = diheral.c;
    //     unsigned int v3 = diheral.d;
    //     std::cout<<"dihedral tag "<<dtag<<" index "<< dihedrali << "reverse index " << dihedral_rtags.data[dtag] <<" vertex:"<<v0<<","<<v1<<","<<v2<<","<<v3<<std::endl;
    // }
}
std::vector< std::string > BondFlipUpdater::getProvidedLogQuantities()
    {
    std::vector<std::string> ret;
    ret.push_back("accepted_bonds");
    m_loggable_quantities = ret; // deep copy. 
    return ret;
    }

Scalar BondFlipUpdater::getLogValue(const std::string& quantity, unsigned int timestep)
    {
    if (std::find(m_loggable_quantities.begin(),m_loggable_quantities.end(), quantity) != m_loggable_quantities.end())
        {
        if (quantity == "accepted_bonds")
            return m_accepted;
        }
    else
        {
        m_exec_conf->msg->error() << "update.bondflip updater: " << quantity
                                  << " is not a valid log quantity." << std::endl;
        throw std::runtime_error("Error getting log value");
        }
    return Scalar(0.0);
    }
std::map<unsigned int, std::vector<unsigned int> > BondFlipUpdater::get_bond_neigh_bonds()
    {
        return m_bond_neigh_bonds;
    }
std::map<unsigned int, std::vector<unsigned int> > BondFlipUpdater::get_bond_faces()
    {
        return m_bond_faces;
    }
std::map<unsigned int, unsigned int> BondFlipUpdater::get_bond_dihedrals()
    {
        return m_bond_dihedral;
    }
/* Export the CPU updater to be visible in the python module
 */
void export_BondFlipUpdater(pybind11::module& m)
    {
    pybind11::class_<BondFlipUpdater, std::shared_ptr<BondFlipUpdater> >(m, "BondFlipUpdater", pybind11::base<Updater>())
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
        std::shared_ptr<ComputeThermo>,
        std::shared_ptr<Logger>,
        std::map<unsigned int, std::vector<unsigned int>>,
        std::map<unsigned int, unsigned int>,
        std::map<unsigned int, std::vector<unsigned int>> >())
        .def("set_params", &BondFlipUpdater::set_params)
        .def("get_bond_neigh_bonds", &BondFlipUpdater::get_bond_neigh_bonds)
        .def("get_bond_dihedrals", &BondFlipUpdater::get_bond_dihedrals)
        .def("get_bond_faces", &BondFlipUpdater::get_bond_faces)
    ;
    }

// ********************************
// here follows the code for BondFlipUpdater on the GPU
// WARNING: This code is not complete! 
#ifdef ENABLE_CUDA

/*! \param sysdef System to perform bond flip
*/
BondFlipUpdaterGPU::BondFlipUpdaterGPU(std::shared_ptr<SystemDefinition> sysdef,
                            std::shared_ptr<ComputeThermo> thermo,
                            std::shared_ptr<Logger> logger,
                            std::map<unsigned int, std::vector<unsigned int>> bond_neigh_bonds,
                            std::map<unsigned int, unsigned int> bond_dihedral,
                            std::map<unsigned int, std::vector<unsigned int>> bond_faces)
        : BondFlipUpdater(sysdef, thermo, logger, bond_neigh_bonds, bond_dihedral, bond_faces)


    {
    m_exec_conf->msg->warning() << "BondFlipUpdaterGPU is not implemented!";
    }


/*! BondFlipUpdaterGPU is not implemented, the original example_plugin content is kept here.
    \param timestep Current time step of the simulation
*/
void BondFlipUpdaterGPU::update(unsigned int timestep)
    {
       if (m_prof) m_prof->push("BondFlipUpdater");

        // access the particle data arrays for writing on the GPU
        ArrayHandle<Scalar4> d_vel(m_pdata->getVelocities(), access_location::device, access_mode::readwrite);

        // call the kernel defined in BondFlipUpdater.cu
        gpu_zero_velocities(d_vel.data, m_pdata->getN());

        // check for error codes from the GPU if error checking is enabled
        if(m_exec_conf->isCUDAErrorCheckingEnabled())
            CHECK_CUDA_ERROR();

        if (m_prof) m_prof->pop();
    }

/* Export the GPU updater to be visible in the python module
 */
void export_BondFlipUpdaterGPU(pybind11::module& m)
    {
    pybind11::class_<BondFlipUpdaterGPU, std::shared_ptr<BondFlipUpdaterGPU> >(m, "BondFlipUpdaterGPU", pybind11::base<BondFlipUpdater>())
        .def(pybind11::init<std::shared_ptr<SystemDefinition>,
        std::shared_ptr<ComputeThermo>,
        std::shared_ptr<Logger>,
        std::map<unsigned int, std::vector<unsigned int>>,
        std::map<unsigned int, unsigned int>,
        std::map<unsigned int, std::vector<unsigned int>> >())
    ;
    }

#endif // ENABLE_CUDA	
