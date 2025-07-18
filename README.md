# pcMMP

Proteome-constrained of *Methanococcus maripaludis* (**pcMMP**). 

## Requirements

Make sure you have **Python 3.11.13+** installed. Required Python packages include:

- cobra==0.29.1  
- matplotlib==3.10.3  
- numpy==2.3.1  
- pandas==2.3.0  

You can install all dependencies using:

  pip install -r requirements.txt

---
## Solver
we used  **SoPlex** as the linear programming (LP) solver (http://soplex.zib.de/).

## Installation

Follow these steps to set up the simulation environment:

#### 1. Clone the Repository

    https://github.com/GhadaSamir2/pcMMP-model.git

#### 2. Navigate to the project directory
    cd pcMMP-model

#### 3. Make the SoPlex binary executable
    chmod +x ./soplex-2.0.0.linux.x86_64.gnu.opt
#### 4. Verify SoPlex is working
    ./soplex-2.0.0.linux.x86_64.gnu.opt

#### 5. Run the notebooks
Open the Jupyter notebooks in your browser to start the simulation


