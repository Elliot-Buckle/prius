from constants import *
from engine import Engine
from nozzle import Nozzle
import numpy as np
from numpy import *
from ocp_vscode import *
import matplotlib.pyplot as plt


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
        
        nozzle.model(export=False, resolution=32)
        
        # Determine masses of cells in array
        
        # Create 2D array for cell masses
        grid_shape = np.shape(nozzle.grid_x)
        print(grid_shape)
        cell_masses = np.full((grid_shape[0] - 1, grid_shape[1] - 1), 'nan')
        
        print(cell_masses)
        
        # Iterate over columns and rows
        for column in range(np.size(cell_masses)):
            for cell in range(column - 1):
                cell_volume = (nozzle.grid_x[cell - 1] + nozzle.grid_x[cell]) * ((nozzle.grid_y[cell] - nozzle.grid_y[cell - 1]) * pi)
                cell = cell_volume / self.material_rho
        
        print("-----------------------")  
        print(cell_masses)
        
if __name__ == "__main__":
    #nozzle = Nozzle(Pc=3*10**6, Tc=3391.91, thrust=300, M=0.026041, mix_ratio=5.3, y=1.2593, nozzle_OD=75*10**-3)
    engine = Engine(fuel="HDPE", oxidizer="N2O", Pc=1.5*10**6, thrust=300, ox_den=786.6, cd=0.44,
                tank_pressure=50.525*10**5, crit_pressure_drop=17.4*10**5, orifice_diameter=1*10**-3, ox_flux=200,
                grain_OD=65*10**-3, cap_OD=75*10**-3, lip_t=10*10**-3, plate_t=10*10**-3, orifice_length = 10*10**-3, manifold_length = 20*10**-3,
                sheath_l=15*10**-3)
    thermal = Thermal(engine.nozzle, k_x=237, k_y=237, material_Cp=1507.248, melting_point=933.5, material_density=2710)
                  
        