#go through a folder containing cif files and extract the residue around a ligand with a cirten cut off
#the ligand HET code must be provided 
import os
import numpy as np
from Bio.PDB import *
from Bio.PDB import PDBIO
from Bio.PDB import Select 

input_dir = input("enter the directory of the .cif files: ")
output_dir = input("enter the directory of the out put files: ")

HET_codes = input("enter list of HET code: ")
cutoff = input("enter the cut off distance: ")


#the location of the "cif" files
os.chdir(input_dir)

file_names = [fn for fn in os.listdir(input_dir) if fn.endswith(".cif")]


for fn in file_names:
    x = 1
    os.chdir(input_dir)
    parser = MMCIFParser()
    structure = parser.get_structure(fn, fn)
    residue = structure[0]

    itrate_res = []
    KKK = None
    for rez in structure[0].get_residues():
       
       if rez.resname in HET_codes and not itrate_res: 
          KKK = rez
          itrate_res.append(KKK)
		  
		  
          binding_residues = []
          residues = list(structure[0].get_residues())
          for res in residues:
            
            if res == KKK:
              continue
            elif res.id[0].startswith("H"):
              continue
            elif KKK == None:
              break
            else:
                try:
                    res_atoms = list(res.get_atoms())
                except KeyError: 
                  continue
                distances = []
                for atom in KKK:
                   for atome in res_atoms:
                        diff_vector = atome.coord - atom.coord
                        distances.append(np.sqrt(np.sum(diff_vector * diff_vector)))
        
        
                for dist in distances:
                     if dist <= cutoff:
                         binding_residues.append(res)
            
          binding_residues.append(KKK)
          class Pocket(Select):
              def accept_residue(self, residue):
                  if residue in binding_residues:
                      return 1
                  else:
                      return 0
          
            
          n = str(rez)
          c = n[9:12]
          fp = fn.replace(".cif", ("_" + str(c) + ".pdb"))      
          os.chdir(output_dir)                
          io = PDBIO()
          io.set_structure(structure)
          
          io.save(fp, Pocket())
          print(str(x) + fp)
          x = x + 1
       else:
           continue 
       

	   
        
    


