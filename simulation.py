#!/usr/bin/env python
# coding: utf-8

from constraints import Lp_File_Constraints
import pandas as pd
import subprocess
import re 
from cobra import*
import cobra
import os
class Simulation:
    
    def __init__(self):
        self.model=None
        self.rxn_blocked=None
        self.lp_file_constraints = Lp_File_Constraints(self.model)
        self.objective=None
        self.optimization_type=None
        self.up_protein_ID=None
        self.up_MW=None
        self.ribosome_ID=None
        self.rRNA_MW=None
        self.mu=None
        self.biomass_ID=None
        self.GAM_ID=None
        self.NGAM_ID=None
        self.NGAM_value=None
        self.kcat_ribo=None
        self.k_mRNA=None
        self.k_deg=None
        self.k_tRNA=None
        self.file_path=None
        self.output_file=None
        self.glycogen_MW=None
        self.glycogen_ID=None
        self.complex_ATP_MW=None
        self.complex_ATP_ID=None
        self.set_constraint=[]
        self.command=None
        
    def read_excel_file(self, filename):
        excel_data = pd.read_excel(filename, sheet_name=None)
        self.lp_file_constraints.complexes_sheet = excel_data['complexes']
        self.lp_file_constraints.Rxn_sheet = excel_data['kcat_values']
        self.lp_file_constraints.gene_protein_sheet = excel_data['gene_protein']
        self.lp_file_constraints.methanogenesis_sheet = excel_data['methanogenesis']
        self.lp_file_constraints.ribosome_sheet = excel_data['ribosome']
        self.lp_file_constraints.tRNA_sheet = excel_data['tRNA']
        self.lp_file_constraints.mRNA_sheet = excel_data['mRNA']
        self.rxn_closed = excel_data['closed_reactions']
     
        
        
    def build_LP_format(self):
        if not os.path.exists("output"):
            os.makedirs("output")
        reaction_ids = [row["rxn_id"] for _, row in self.rxn_closed.iterrows()]
        self.close_reactions(reaction_ids)
        self.lp_file_constraints.model=self.model
        self.lp_file_constraints.mu=self.mu
        self.lp_file_constraints.open_lp_file(self.file_path,self.objective,self.optimization_type)
        self.lp_file_constraints.apply_sv_constraint()
        self.lp_file_constraints.apply_kcat_constrain()
        self.lp_file_constraints.metabolic_mass()
        self.lp_file_constraints.up_mass(self.up_protein_ID, self.up_MW)
        self.lp_file_constraints.apply_ribosome_constrain(self.ribosome_ID, self.kcat_ribo)
        self.lp_file_constraints.rRNA_mass(self.ribosome_ID, self.rRNA_MW)
        self.lp_file_constraints.ribosome_mass(self.ribosome_ID)
        self.lp_file_constraints.methanogenesis_mass( )
        self.lp_file_constraints.glycogen_mass(self.glycogen_MW,self.glycogen_ID)
        self.lp_file_constraints.ATP_mass(self.complex_ATP_MW,self.complex_ATP_ID)
        self.lp_file_constraints.apply_tRNA_constraint(self.k_tRNA)
        self.lp_file_constraints.tRNA_mass()
        self.lp_file_constraints.apply_mRNA_constraint(self.k_mRNA,self.k_deg)
        self.lp_file_constraints.mRNA_mass(self.k_deg)
        self.lp_file_constraints.add_constraint("metaboicMass + ribosome_mass + UP - total_protein = 0")
        
        self.lp_file_constraints.add_constraint("rRNA + tRNA + mRNA  - RNA= 0")
       
        self.lp_file_constraints.add_constraint("total_protein + RNA + Glycogen = 0.93")
        
        if self.set_constraint:
           
            for constraint in self.set_constraint:
                
                self.lp_file_constraints.add_constraint(constraint)
            
        self.lp_file_constraints.add_bounds( self.biomass_ID,self.GAM_ID,self.NGAM_ID,self.NGAM_value)
        self.lp_file_constraints.close_lp_file()
        
        
    def close_reactions(self,reaction_ids):
        for reaction_id in reaction_ids:
            reaction=self.model.reactions.get_by_id(reaction_id)
            reaction.lower_bound=0
            reaction.upper_bound=0
   
    def GAM(self,GAM_ID,GAM):
        
        reaction=self.model.reactions.get_by_id(GAM_ID)
        reaction.reaction=f" {GAM} ATP[c0] + {GAM} H2O[c0] ==>  {GAM} ADP[c0] + {GAM} Phosphate[c0] + {GAM} H[c0]"
        

    def check_mu_optimal(self,mu):
        self.mu=mu
        self.build_LP_format() 
        command=f"./soplex-2.0.0.linux.x86_64.gnu.opt -s0 -x -q -c -f1e-20 -o1e-20 --solvemode=2 --readmode=1  --writebas=basefile.bas {self.file_path} > {self.output_file}"
        
        with open(self.output_file, "w") as out:
            subprocess.run(command, shell=True, stdout=out, stderr=subprocess.DEVNULL)
        
        output_file = open(self.output_file,"r")
        output_file=output_file.read()
        objective_value = re.search(r"SoPlex status\s+:\s+problem is solved \[([^\]]+)\]", output_file)
        if objective_value:
            objective_value = objective_value.group(1)
            if objective_value == "optimal":
                return True
            else:
                return False
            
    def search_max_growth_rate(self,start, end):
        increment = 0.001
        max_growth_rate = None
        while start <= end:
            mid = (start + end) / 2
            if self.check_mu_optimal(mid):
                max_growth_rate = mid
                start = mid + increment
            else:
                end = mid - increment

        if max_growth_rate:
            self.check_mu_optimal(max_growth_rate)
        return max_growth_rate


    def set_parameter(self,GAM=None,NGAM=None, objective=None, optimization_type="Minimize",
                      biomass_ID="R4799", GAM_ID="R4800", NGAM_ID="R4801", glycogen_ID="R2879", glycogen_MW=162,
                      up_protein_ID="R2871", up_MW=49919.951183051555, ribosome_ID="R2061", rRNA_MW=1444610,k_mRNA=0.034,k_deg=8.3,
                      k_tRNA=2.73,kcat_ribo=22, file_path=None, output_file=None, model_file=None, 
                      excel_file=None):
        
        self.model = cobra.io.read_sbml_model(model_file)
        self.objective = objective
        self.optimization_type = optimization_type
        self.biomass_ID = biomass_ID
        self.GAM_ID = GAM_ID
        self.NGAM_ID = NGAM_ID
        self.NGAM_value=NGAM
        self.glycogen_ID = glycogen_ID
        self.glycogen_MW = glycogen_MW
        self.complex_ATP_ID = "R2730"
        self.complex_ATP_MW = 1265596.72236614
        self.up_protein_ID = up_protein_ID
        self.up_MW = up_MW
        self.ribosome_ID = ribosome_ID
        self.rRNA_MW = rRNA_MW
        self.kcat_ribo = kcat_ribo
        self.k_mRNA=k_mRNA
        self.k_deg=k_deg
        self.k_tRNA=k_tRNA
        self.file_path = file_path
        self.output_file = output_file
        self.read_excel_file(excel_file)
        if GAM is not None :
            self.GAM(GAM_ID, GAM)