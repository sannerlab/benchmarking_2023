#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 13:29:31 2022
v2
@author: sshanker

Information:
Collection of functions for the OpenMM based restrained minimization,
and interaction energy calculation.   

"""

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from pdbfixer import PDBFixer
import prody
from MolKit2 import Read
from prody.measure.contacts import findNeighbors


def split_pdb_to_chain_A_and_B(pdbid):    
    '''This prorams expect chain ID for receptor as 'A' 
    and for peptide as 'B'
    '''
    out_pdb_int = pdbid[:-4]    
    mol = Read(pdbid)
    
    rec = mol._ag.select('chid A')
    pep = mol._ag.select('chid B')
    
    prody.writePDB(out_pdb_int+"_A.pdb", rec)
    prody.writePDB(out_pdb_int+"_B.pdb", pep)


def identify_interface_residues(pdb_file, needs_b=0):
    '''This prorams expect chain ID for receptor as 'A' 
    and for peptide as 'B'
    '''
    mol_temp = Read(pdb_file)
    
    rec = mol_temp._ag.select('chid A not hydrogen')
    pep = mol_temp._ag.select('chid B not hydrogen')
    

    near_n = findNeighbors(rec, 5 , pep)

    interacting_chainA_res = []
    if needs_b ==1:
        interacting_chainB_res = []
    for a1,a2,d in near_n:
        
        interacting_chainA_res.append(a1.getResindex())
        if needs_b ==1:
            interacting_chainB_res.append(a2.getResindex())        
        
    interacting_chainA_res = list(set(interacting_chainA_res))
    
    interacting_chainA_res.sort()
    # import pdb; pdb.set_trace()
    res_list_prody=[]
    for i in interacting_chainA_res:
        res_list_prody.append(rec.select(' name CA and resindex %d' % i ).getResnames()[0])
    print("prody ignore:",res_list_prody)
    
    if needs_b ==1:
        interacting_chainB_res = list(set(interacting_chainB_res))
        interacting_chainB_res.sort()
        pep_list_prody=[]
        for i in interacting_chainB_res:
            pep_list_prody.append(pep.select(' name CA and resindex %d' % i ).getResnames()[0])
        print("prody ignore(B):",pep_list_prody)
        
        return interacting_chainA_res,interacting_chainB_res
   
    # interacting_chainA_res=np.array(interacting_chainA_res)
    return interacting_chainA_res
    


def fix_my_pdb(pdb_in,out=None):
    '''This program fixes pdbs for simple errors
    like wrong atom name assignment from AD. 
    It also uses PDBFixer to add missing atoms'''
    if out==None:
        pdb_out = "./fixed/" + pdb_in.split('/')[-1].split('.')[0]+'_fixed.pdb'
    else:
        pdb_out = out
        
    mol_tmp =Read(pdb_in)

    rs_names = mol_tmp._ag.getResnames()
    atom_names = mol_tmp._ag.getNames()
    element_names = mol_tmp._ag.getElements()

    for indx, (i,j) in enumerate(zip(rs_names, atom_names)):
        if i == 'LEU':
            if j == 'CD':
                atom_names[indx]='CD1'
                element_names[indx]='C'
            elif j == 'CG':
                element_names[indx]='C'
                
        elif i == 'ASN':
            if j == '2HD':
                atom_names[indx]='2HD2'
                element_names[indx]='H'
        elif i == 'ARG':
            if j == '2HN1':
                atom_names[indx]='2HH1'
                element_names[indx]='H'
        elif i == 'VAL':
            if j == 'CD1':
                atom_names[indx]='CG1'
                element_names[indx]='C'
            elif j == 'CD2':
                atom_names[indx]='CG2'
                element_names[indx]='C'
       
                
            
    mol_tmp._ag.setNames(atom_names)
    mol_tmp._ag.setElements(element_names)
    prody.writePDB(pdb_out,mol_tmp._ag)
        
    fixer = PDBFixer(filename=pdb_out)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.findNonstandardResidues()
    
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    fixer.removeHeterogens(False)
    
    with open( pdb_out , 'w') as outfile:
         PDBFile.writeFile(fixer.topology, fixer.positions, file=outfile,
    keepIds=True)
    
    return pdb_out

  
def restrain_(system, pdb, chain='A',ignore_list=[]):
    '''To apply strong harmonic restrains on non-interface residues'''
    restraint = CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")    
    restraint.addGlobalParameter('k', 20000.0*kilojoules_per_mole/nanometer*nanometer)
    
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')
    # print(ignore_list)
    res_list_openmm=[]
    for atom in pdb.topology.atoms():
        if len(ignore_list) > 0:
            if ignore_list.count(atom.residue.index) > 0:
                if atom.name == 'CA':
                    res_list_openmm.append(atom.residue.name)
                continue
        
        # import pdb ; pdb.set_trace()
        if chain =='A':
            if atom.element.name == 'hydrogen':
                continue
            if atom.residue.chain.id == 'A':       
                restraint.addParticle(atom.index, pdb.positions[atom.index])
        else:
            restraint.addParticle(atom.index, pdb.positions[atom.index])
    if chain =='A':       
        print("omm ignore  :",res_list_openmm)            
    system.addForce(restraint)
    
    
    
def get_single_point_energy(pdbfile,mode='vacuum'):
    'To get a single point energy of a pdb file'
    pdb_handle = PDBFile(pdbfile)

    if mode == 'implicit':
        force_field = ForceField("amber99sb.xml",'implicit/gbn2.xml')
        print("Using GBSA (gbn2) environment for energy calculation")
    
    else:
        force_field = ForceField("amber99sb.xml")
        print("Using In-Vacuo environment for energy calculation")
        
    system = force_field.createSystem(pdb_handle.topology, nonbondedCutoff=1*nanometer, constraints=HBonds)
    restrain_(system, pdb_handle, 'All')
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    simulation = Simulation(pdb_handle.topology, system, integrator)
    simulation.context.setPositions(pdb_handle.positions)
    simulation.minimizeEnergy(maxIterations=1)
    state = simulation.context.getState(getEnergy=True)
    # import pdb; pdb.set_trace()
        
    return state.getPotentialEnergy()._value

          

def _openmm_minimize( pdb_str: str,env='implicit'):
    """Minimize energy via openmm.
    Adopted from AF2 minimization code """    
    
    pdb = PDBFile(pdb_str)
    
    if env == 'implicit':
        force_field = ForceField("amber99sb.xml",'implicit/gbn2.xml')
        print("Using GBSA (gbn2) environment for energy minimization")
    else:
        force_field = ForceField("amber99sb.xml")
        print("Using in-vacuo environment for energy minimization")
    # integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    constraints = HBonds
    
    system = force_field.createSystem(  pdb.topology, nonbondedCutoff=1*nanometer, 
                                      constraints=constraints)
    
    
    ignore_list = identify_interface_residues(pdb_str) # No restrained residues
    restrain_(system, pdb, ignore_list=ignore_list)


    # platform = openmm.Platform.getPlatformByName("CUDA" if use_gpu else "CPU")
    simulation = Simulation( pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
      
    ret = {}
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    ret["einit"] = state.getPotentialEnergy()
    # ret["posinit"] = state.getPositions(asNumpy=True)
    simulation.minimizeEnergy(maxIterations=100, tolerance=0.01)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    ret["efinal"] = state.getPotentialEnergy()
    # ret["pos"] = state.getPositions(asNumpy=True)
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(pdb_str[:-4]+"_min.pdb", 'w'))
   
    return ret,system



def estimate_energies_for_pdb(pdb_file, env='implicit'):
    print('working on:',pdb_file.split("/")[-1])
    fixed_pdb = fix_my_pdb(pdb_file, "./fixed_general/"+ pdb_file.split("/")[-1][:-4] +"_fixed.pdb")
    print("minimizing ...")
    _=_openmm_minimize(fixed_pdb,env)
    minimized_pdb = fixed_pdb[:-4]+"_min.pdb"
    enzs=[] # energy values
    split_pdb_to_chain_A_and_B(fixed_pdb[:-4]+"_min.pdb")
    for pdbs in [minimized_pdb, minimized_pdb[:-4]+"_A.pdb",minimized_pdb[:-4]+"_B.pdb" ]:
        enzs.append(get_single_point_energy(pdbs,env))

    enzs.append(enzs[0] - enzs[1] -enzs[2])
    enzs.append(enzs[0] - enzs[1])
    e_str=''
    for i in enzs:
        e_str = e_str + "%10.2f" % i
    print(' E_Complex E_Receptr E_Peptide dE_Interaction dE_ComRes')   
    print (e_str)
    return enzs
    

