# @Author: Siyu Li(sli032@ucr.edu)
# @Date:   2022-06-07 15:39:41
# @Last Modified by:   Siyu Li
# @Last Modified time: 2022-06-07 15:42:40
import numpy as np
import pandas as pd
import json

ONEDGE=32767


class mem():
    """docstring for Shell"""
    def __init__(self):
        # super(polymer, self).__init__()
        self.vertex = pd.DataFrame(columns = ['x','y','z','l','t','edge'])
        self.line = pd.DataFrame(columns = ['v0','v1','t','edge']).astype('int32')
        self.triangle = pd.DataFrame(columns = ['v0','v1','v2','l0','l1','l2','edge']).astype('int32')
        self.n = 0
        self.dlist = []
        self.ptypelist = []
        self.particles = []
        self.diameter = 1.0
        self.bodyid = -1
        self.defaulttype = 'M'
        self.btypename = 'mbond'
        self.dtypename = 'mdihedral'


    def other_vertex(self,t,v0,v1):
        [tv0, tv1, tv2] = self.triangle.loc[t, ['v0','v1','v2']]
        return tv0 + tv1 + tv2 - v0 - v1

    # def update_particles(self):
    #     self.n = len(self.vertex)
    #     self.dlist = [self.diameter] * self.n
    #     self.ptypelist = [self.defaulttype] * self.n
    #     self.particles = np.array(self.vertex[list('xyz')])

    # def update_bonds(self):
    #     self.bonds = np.int_(self.line[self.line.edge > 0][['v0','v1']])
    
    # def update_angles(self):
    #     self.angles = np.int_(self.triangle[['v0', 'v1', 'v2']])
        
    # def update_dihedrals(self):
    #     dihedralgroup = []
    #     for l in self.line[self.line.edge!=ONEDGE].index:
    #         l_pair = self.line.loc[l,'edge']
    #         if l_pair > 0:
    #             [v0, v1, t0] = self.line.loc[l, ['v0','v1','t']]
    #             [v10, v11, t1] = self.line.loc[l_pair, ['v0','v1','t']]
    #             v2 = self.other_vertex(t0, v0, v1)
    #             v12 = self.other_vertex(t1, v10, v11)
    #             dihedralgroup = dihedralgroup + [[v2, v1, v11, v12]]
    #     self.dihedrals = np.int_(dihedralgroup)

    # def update_mem_info(self):
    #     self.update_particles()
    #     self.update_bonds()
    #     self.update_angles()
    #     self.update_dihedrals()
    #     print(self.particles)
    #     print(self.bonds)
    #     print(self.angles)
    #     print(self.dihedrals)
    
    # def update_bond_flip_info(self):
    #     bond_neigh_bonds, bond_dihedrals, bond_faces = {}, {}, {}
    #     bond_dic, tri_dic = {}, {}
        
    #     for i, l in enumerate(self.line[self.line.edge>0].index):
    #         bond_dic[l] = i
    #     for i, t in enumerate(self.triangle.index):
    #         tri_dic[t] = i
    #     dihedral_index = 0
    #     face_index = 0
    #     for l in self.line[self.line.edge!=ONEDGE].index:
    #         l_pair = self.line.loc[l,'edge']
    #         if l_pair > 0:
    #             t0 = self.line.loc[l,'t']
    #             t1 = self.line.loc[l_pair,'t']
    #             lines0 = set(self.line[self.line.t==t0].index)
    #             lines1 = set(self.line[self.line.t==t1].index)
    #             neigh_bonds = list(lines0.union(lines1) - set([l, l_pair]))
    #             bond_index = bond_dic[l]
    #             bond_neigh_bonds[bond_index] = []
    #             for bond in neigh_bonds:
    #                 bondi = min(bond, abs(self.line.loc[bond, 'edge']))
    #                 bond_neigh_bonds[bond_index].append(bond_dic[bondi])
    #             bond_dihedrals[bond_index] = dihedral_index
    #             dihedral_index += 1
    #             face_index += 2
    #     self.bond_neigh_bonds = bond_neigh_bonds
    #     self.bond_dihedrals = bond_dihedrals
        
    #     self.bond_faces={}
    #     for l in self.line[self.line.edge!=ONEDGE].index:
    #         l_pair = self.line.loc[l,'edge']
    #         if l_pair > 0:
    #             [v0, v1] = self.line.loc[l, ['v0', 'v1']]
    #             t0, t1 = self.line.loc[l, 't'], self.line.loc[l_pair, 't']
    #             self.bond_faces[bond_dic[l]] = [tri_dic[t0], tri_dic[t1]]
    #     print("bond_faces:", self.bond_faces)
    
        
    # def print(self):
    #     cprint('Shell Structure:\n','yellow')
    #     print(self.vertex)
    #     print(self.line)


    def update_particles(self):
        self.n=len(self.vertex)
        self.dlist=[self.diameter]*self.n
        self.ptypelist = [self.defaulttype] * self.n
        self.pdic=dict(zip(self.vertex.index, range(self.n)))
        self.particles=np.array(self.vertex[list('xyz')])

    def correct_index(self,obj):
        for i,row in enumerate(obj):
            obj[i]=list(map(lambda x: self.pdic[x], row))
        return np.array(obj)

    def update_bonds(self):
        bonds = np.int_(self.line[self.line.edge > 0][['v0','v1']])
        self.bonds=self.correct_index(bonds)        
    
    def update_angles(self):
        angles = np.int_(self.triangle[['v0', 'v1', 'v2']])
        self.angles=self.correct_index(angles) 
        
    def update_dihedrals(self):
        dihedralgroup=[]
        for l in self.line[self.line.edge!=ONEDGE].index:
            l_pair=self.line.loc[l,'edge']
            if l_pair>0:
                [v0,v1,t0]=self.line.loc[l,['v0','v1','t']]
                [v10,v11,t1]=self.line.loc[l_pair,['v0','v1','t']]
                v2=self.other_vertex(t0,v0,v1)
                v12=self.other_vertex(t1,v10,v11)
                dihedralgroup=dihedralgroup+[[v2,v1,v11,v12]]
        self.dihedrals=self.correct_index(dihedralgroup)

    def update_mem_info(self):
        self.update_particles()
        self.update_bonds()
        self.update_angles()
        self.update_dihedrals()
        print(self.particles)
        print(self.bonds)
        print(self.angles)
        print(self.dihedrals)
    
    def update_bond_flip_info(self):
        bond_neigh_bonds, bond_dihedrals, bond_faces = {}, {}, {}
        bond_dic, tri_dic = {}, {}
        
        for i, l in enumerate(self.line[self.line.edge>0].index):
            bond_dic[l] = i
        for i, t in enumerate(self.triangle.index):
            tri_dic[t] = i
        dihedral_index = 0
        face_index = 0
        for l in self.line[self.line.edge!=ONEDGE].index:
            l_pair=self.line.loc[l,'edge']
            if l_pair>0:
                t0 = self.line.loc[l,'t']
                t1 = self.line.loc[l_pair,'t']
                lines0 = set(self.line[self.line.t==t0].index)
                lines1 = set(self.line[self.line.t==t1].index)
                neigh_bonds = list(lines0.union(lines1) - set([l, l_pair]))
                bond_index = bond_dic[l]
                bond_neigh_bonds[bond_index] = []
                for bond in neigh_bonds:
                    bondi = min(bond, abs(self.line.loc[bond, 'edge']))
                    bond_neigh_bonds[bond_index].append(bond_dic[bondi])
                bond_dihedrals[bond_index] = dihedral_index
                dihedral_index += 1
                face_index += 2
        self.bond_neigh_bonds = bond_neigh_bonds
        self.bond_dihedrals = bond_dihedrals
        
        self.bond_faces={}
        for l in self.line[self.line.edge!=ONEDGE].index:
            l_pair=self.line.loc[l,'edge']
            if l_pair>0:
                [v0, v1] = self.line.loc[l, ['v0', 'v1']]
                t0, t1 = self.line.loc[l, 't'], self.line.loc[l_pair, 't']
                self.bond_faces[bond_dic[l]] = [tri_dic[t0], tri_dic[t1]]
        print("bond_faces:", self.bond_faces)
        
    def print(self):
        cprint('Shell Structure:\n','yellow')
        print(self.vertex)
        print(self.line)

    def frame_in(self, file):
        with open('./{}.csv'.format(file),'r') as f:
            frame_row = 0
            for row_num, line in enumerate(f, 1):
                if line[0] == '#':
                    data_row = row_num
                    data_name, data_nrow = line.split(':')
                    print(data_name[1:],int(data_nrow))
                    df = pd.read_csv('./{}.csv'.format(file), skiprows = data_row, index_col = 0,header = 0, nrows = int(data_nrow),error_bad_lines = False)
                    print(len(df))
                    print(df)
                    if data_name[1:] == 'vertex':
                        self.vertex = df
                    elif data_name[1:] == 'line':
                        self.line = df
                    elif data_name[1:] == 'triangle':
                        self.triangle = df
                    else: raise NameError('Unknown input data')