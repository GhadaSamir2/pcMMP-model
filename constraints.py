import itertools
from ast import literal_eval
import pandas as pd
import cobra
from cobra import*

class Lp_File_Constraints:
    def __init__(self,model):
        self.excel_data = None
        self.complexes_sheet = None
        self.Rxn_sheet = None
        self.gene_protein_sheet = None
        self.methanogenesis_sheet = None
        self.ribosome_sheet = None
        self.tRNA_sheet = None
        self.mRNA_sheet = None
        self.model=model
        self.lp_file=None
        self.mu=None
        self.counter = None

        
    def open_lp_file(self,file_path,objective,optimization_type):
        
        self.counter=itertools.count(start=0)
        self.lp_file = open(file_path, "w")
        self.lp_file.write("{} \n".format(optimization_type))
        self.lp_file.write("obj: {} \n".format(objective))
        self.lp_file.write("Subject To\n")
        
       

    def apply_sv_constraint(self):
        s_matrix = cobra.util.array.create_stoichiometric_matrix(self.model)
        num_metabolites, num_reactions = s_matrix.shape
        for i in range(num_metabolites):
            equation = ""
            c = next(self.counter)
            equation = f'c{c}: '
            for j in range(num_reactions):
                if s_matrix[i][j] != 0:
                    if s_matrix[i][j] > 0:
                        equation += " + "
                    elif s_matrix[i][j] < 0:
                        equation += " - "
                    equation += f'{abs(s_matrix[i][j])} R{j + 1}'
                if len(equation) >= 200:
                    self.lp_file.write(equation + '\n')
                    equation = ""
            equation = f'{equation} = 0\n'
            self.lp_file.write(equation)

    def apply_kcat_constrain(self):
        for _, row in self.complexes_sheet.iterrows():
            Reactions = literal_eval(row["Reactions"])
            if len(Reactions) > 1:
                complex_formation_ID = row["complex_ID"]
                c = next(self.counter)
                self.lp_file.write("c{}: ".format(c))
                for reaction_ID in Reactions:
                    kcat = float(self.Rxn_sheet[self.Rxn_sheet['rxn_id'] == reaction_ID].kcat.iloc[0])
                    if reaction_ID!="R1161":
                        if kcat < 0.75:
                            kcat = 0.75
                    kcat = (1 / (3600 * kcat))
                    self.lp_file.write(" + {} {}".format(kcat, reaction_ID))
                M = (1 / self.mu)
                self.lp_file.write(" - {} {} = 0\n".format(M, complex_formation_ID))
            else:
                reaction = literal_eval(row["Reactions"])
                reaction_ID = reaction[0]
                complex_formation_ID = row["complex_ID"]
                kcat = float(self.Rxn_sheet[self.Rxn_sheet['rxn_id'] == reaction_ID].kcat.iloc[0])
                if kcat < 0.75:
                    kcat = 0.75
                c = next(self.counter)
                M = (3600 * kcat) * (1 / self.mu)
                self.lp_file.write("c{}: {} - {} {} = 0\n".format(c, reaction_ID, M, complex_formation_ID))

                
                
    def metabolic_mass(self):
        c = next(self.counter)
        equation_protein=f'c{c}: '

        for _,row in self.complexes_sheet.iterrows() :
            complex_dilution_ID=row["ID_dilution"]
            MW=float(row["MW"])
            M=(1/self.mu)*(MW/1000)
            equation_protein+='+ {}  {} '.format(M,complex_dilution_ID)

            if len(equation_protein) >= 200:

                self.lp_file.write(equation_protein + '\n')
                equation_protein = ""

        equation_protein+=' - metaboicMass  = 0\n'.format()   
        self.lp_file.write(equation_protein)

            
    def up_mass(self, up_protein_ID, up_MW):
        c= next(self.counter)
        equation_protein=f'c{c}: '
        M=( 1 / self.mu )
        MW=(up_MW/1000)
        value=M*MW
        equation_protein+='{} {} - UP = 0\n'.format(value,up_protein_ID)   
        self.lp_file.write(equation_protein)

        
        
    def apply_ribosome_constrain(self, ID_ribosome, kcat):
        equation_ribosome=''
        counter_eR=itertools.count(start=1)
        for _,row in self.gene_protein_sheet.iterrows():
            equation=''
            translation_reaction_ID=row['translation_reaction_ID']
            lenght=row['length_peptide']
            kcat_ribo=( ( kcat * 3600 ) * self.mu)
            kcat_per_protein=kcat_ribo/lenght
            c= next(self.counter)
            i=next(counter_eR)

            self.lp_file.write('c{}: {} - {} eR{} = 0 \n'.format(c,translation_reaction_ID,kcat_per_protein,i))

            equation_ribosome+="eR{} + ".format(i)

        c= next(self.counter)
        M=( 1 / self.mu )
        equation_ribosome=equation_ribosome[:-2]
        equation_ribosome+="- {} {}".format(M,ID_ribosome)
        self.lp_file.write("c{}: {} = 0 \n".format(c,equation_ribosome))
    
    

    def rRNA_mass(self, ID_ribosome, rRNA_MW):
        equation_rRNA=''
    
        M=( 1 / self.mu ) * ( rRNA_MW / 1000 )
        c= next(self.counter)
        equation_rRNA=" {} {} - rRNA".format(M,ID_ribosome)

        self.lp_file.write("c{}: {} = 0 \n".format(c,equation_rRNA))
        
        

    def ribosome_mass(self, ID_ribosome):
        
        c= next(self.counter)
        MW=0
        equation_protein=f'c{c}: '
        for _,row in self.ribosome_sheet.iterrows() :

            MW+=float(row["MW"])

        M=( MW / 1000 )
        value=( 1 / self.mu ) * ( M )

        equation_protein+=' {} {} - ribosome_mass = 0\n'.format(value,ID_ribosome)   
        self.lp_file.write(equation_protein)
        

    def methanogenesis_mass(self):
        c= next(self.counter)
        equation_protein=f'c{c}: '
        for _,row in self.methanogenesis_sheet.iterrows() :

            complex_dilution_ID=row["ID_dilution"]
            MW=float(row["MW"])
            M=( MW / 1000 )
            value=( 1 / self.mu ) * ( M )
            equation_protein+='+ {}  {} '.format(value,complex_dilution_ID)

        equation_protein+=' - Methanogenesis_mass = 0\n'.format()   
        self.lp_file.write(equation_protein)
       
        
    def glycogen_mass(self,glycogen_MW,glycogen_ID):
    
        c= next(self.counter)
        equation_glycogen=f'c{c}: '
        M=(1/self.mu)
        MW=(glycogen_MW/1000)
        value=M*MW
        equation_glycogen+='{} {} - Glycogen = 0\n'.format(value,glycogen_ID)   
        self.lp_file.write(equation_glycogen)
    
        
   
    def ATP_mass(self,complex_ATP_MW,complex_ATP_ID):
    
        c= next(self.counter)
        equation=f'c{c}: '
        M=(1/self.mu)
        MW=(complex_ATP_MW/1000)
        value=M*MW
        equation+='{} {} - ATP_Synthase = 0\n'.format(value,complex_ATP_ID)   
        self.lp_file.write(equation)
        
        
    def apply_tRNA_constraint(self,k_tRNA):
       
        for _,row in self.tRNA_sheet.iterrows() :
            c= next(self.counter)
            equation=f'c{c}: '
            ID_charging=row["ID_charging"]
            ID_transcription=row["ID_transcription"]
            value=( 1 / self.mu ) * (k_tRNA*3600 )
            equation='c{} : {} - {}  {} = 0  \n'.format(c,ID_charging,value,ID_transcription)
            self.lp_file.write(equation)
        
        
    def tRNA_mass(self):
        c= next(self.counter)
        equation=f'c{c}: '
        for _,row in self.tRNA_sheet.iterrows() :
            ID=row["ID_transcription"]
            MW=float(row["MW"])
            M=( MW / 1000 )
            value=( 1 / self.mu ) * ( M )
            equation+='+ {}  {} '.format(value,ID)

        equation+=' - tRNA = 0\n'.format()   
        self.lp_file.write(equation)
        
    def apply_mRNA_constraint(self,k_mRNA,k_deg):
       
        for _,row in self.mRNA_sheet.iterrows() :
            c= next(self.counter)
            equation=f'c{c}: '
            ID_translation=row["ID_translation"]
            ID_transcription=row["ID_transcription"]
            value=(k_mRNA*3600)/(k_deg+self.mu )
            equation='c{} : {} - {}  {} = 0  \n'.format(c,ID_translation,value,ID_transcription)
            self.lp_file.write(equation)
            
            
    def mRNA_mass(self,k_deg):
        c= next(self.counter)
        equation=f'c{c} : '
        for _,row in self.mRNA_sheet.iterrows():
            
            transcription_ID=row["ID_transcription"]
           
            MW=float(row["MW"])
            M=( MW / 1000 )
            value=( 1 / (self.mu+k_deg) ) * ( M )
            equation+='+ {}  {} '.format(value,transcription_ID)
            
            if len(equation) >= 200:

                self.lp_file.write(equation + '\n')
                equation= ""

        equation+=' - mRNA = 0\n'.format()   
        self.lp_file.write(equation)
        
      
    def add_constraint(self,constraint):
        c= next(self.counter)
        self.lp_file.write("c{}: {} \n".format(c,constraint))
        

    def add_bounds(self, biomass_ID,GAM_ID,NGAM_ID,NGAM):
        self.lp_file.write("Bounds\n")
        num_reactions=len(self.model.reactions)

        for reaction_index in range(1,num_reactions+1):
            reaction_id=f'R{reaction_index}'

            if reaction_id == biomass_ID:
                self.lp_file.write(" {} <= R{} <= {}  \n".format(self.mu,reaction_index,self.mu))
            elif reaction_id == GAM_ID:
                self.lp_file.write(" {} <= R{} <= {}  \n".format(self.mu,reaction_index,self.mu))
            elif reaction_id==NGAM_ID:
               
                self.lp_file.write(" {} <= R{} <= {}  \n".format(NGAM,reaction_index,NGAM))
                
            else:

                if self.model.reactions.get_by_id(reaction_id).upper_bound==0:
                    
                    self.lp_file.write("0 <= R{} <= 0 \n".format(reaction_index))
                else:   

                    self.lp_file.write("0 <= R{} <= +Inf \n".format(reaction_index))

    
     
    def close_lp_file(self):
        self.lp_file.write("End\n")
        self.lp_file.close()
            
