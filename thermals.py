from constants import *
from nozzle import Nozzle
import numpy as np
from numpy import *
import matplotlib.pyplot as plt

# Thermal simulation
class Thermal:
    def __init__(
        self,
        nozzle:Nozzle,
        k_x:float,
        k_y:float,
        material_Cp:float,
        melting_point:float,
        material_density:float,
    ):
        self.conductivity_x = k_x
        self.conductivity_y = k_y
        self.material_Cp = material_Cp
        self.melting_point = melting_point
        self.material_rho = material_density
        
        # Determine masses of cells in array
        
        # Create 2D array for cell masses
        cell_masses = np.full((nozzle.model_xpoints - 1, nozzle.model_ypoints), 'nan')
        
        # Iterate over columns and rows
        for column in cell_masses:
            for cell in column:
                cell_volume = (nozzle.model_xpoints[cell - 1] + nozzle.model_xpoints[cell]) * ((nozzle.model_ypoints[cell] - nozzle.model_ypoints[cell - 1]) * pi)
                cell = cell_volume / material_density
                
        print(cell_masses)
        
if __name__ == "__main__":
    nozzle = Nozzle(Pc=3*10**6, Tc=3391.91, thrust=300, M=0.026041, mix_ratio=5.3, y=1.2593, nozzle_OD=1)
    thermal = Thermal(nozzle, k_x=237, k_y=237, material_Cp=1507.248, melting_point=933.5, material_density=2710)
                  
        