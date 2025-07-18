#!/usr/bin/env python
# coding: utf-8

import re
import cobra
import itertools
import pandas as pd
class Fluxes:
    def __init__(self,model,output_file):
        self.model=model
        self.output_file=output_file
        self.complex_sheet=None
    
    def extract_reaction_flux(self, reaction_id):
        output_file = open(self.output_file, "r")
        output_file = output_file.read()
        pattern = r"(?<!e){}\s+(\d+(?:\.\d+)?)".format(reaction_id)
        reaction_flux = re.search(pattern, output_file)
        if reaction_flux:
            reaction_flux = reaction_flux.group(1)
            reaction_flux = float(reaction_flux)
        else:
            reaction_flux = 0
        return reaction_flux

    def calculate_v_complex(self,mu):
        num_reactions=len(self.model.reactions)

        reaction_ids=[]
        concentration=[]
        v_proteins=[]
        reactions=[]
        output_file = open(self.output_file,"r")
        output_file=output_file.read()
        for _,row in self.complex_sheet.iterrows():

            reaction_id=row["ID_dilution"]
            reaction_ids.append(reaction_id)
           
            reactions.append(row["Reactions"])

            MW=row["MW"]

            pattern =r"(?<!e){}\s+(\d+(?:\.\d+)?)".format(reaction_id)
            match = re.search(pattern, output_file)
            if match:
                value = match.group(1)
                value = float(value)
                v_per_protein =(value/mu)*(MW/1000)
                v_proteins.append(v_per_protein)
            else:

                v_proteins.append(0)
          
        
        total_protein= self.extract_reaction_flux("total_protein")
        concentration=[x/total_protein for x in v_proteins]
        reactions_data = pd.DataFrame({"complex_ID":reaction_ids,"concentration":concentration,"reactions":reactions})
        return reactions_data
        
    def extract_fluxes(self):
        num_reactions=len(self.model.reactions)
        reactions_data=[]
        for reaction_index in range(1,num_reactions+1):
            reaction_id=f'R{reaction_index}'
            reaction=self.model.reactions.get_by_id(reaction_id)
            reaction_flux=self.extract_reaction_flux(reaction_id)
            reactions_data.append({"reaction_id":reaction.id,"reaction":reaction.reaction,"flux":reaction_flux})
        return reactions_data
   

