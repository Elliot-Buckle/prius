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
        array_x_points = np.ndarray([nozzle.grid_x_points])
        array_y_points = np.ndarray([nozzle.grid_y_points])
        
        # Determine masses of cells in array
        
        # Create 2D array for cell masses
        cell_masses = np.full((array_x_points - 1, array_y_points - 1), 'nan')
        
        print(cell_masses)
        
        # Iterate over columns and rows
        for column in cell_masses:
            for cell in column:
                cell_volume = (array_x_points[cell - 1] + array_x_points[cell]) * ((array_y_points[cell] - array_y_points[cell - 1]) * pi)
                cell = cell_volume / self.material_rho
                
        print(cell_masses)
        
if __name__ == "__main__":
    #nozzle = Nozzle(Pc=3*10**6, Tc=3391.91, thrust=300, M=0.026041, mix_ratio=5.3, y=1.2593, nozzle_OD=75*10**-3)
    engine = Engine(Pc=1.5*10**6, Tc=3252.15, thrust=300, M=0.025423, mix_ratio=6.1, y=1.2501, ox_den=786.6, cd=0.44, viscosity=3.237*10**-4,
                tank_pressure=50.525*10**5, crit_pressure_drop=17.4*10**5, orifice_diameter=1*10**-3, a=0.417, n=0.347, fuel_den=902, ox_flux=200,
                grain_OD=65*10**-3, cap_OD=75*10**-3, lip_t=10*10**-3, plate_t=10*10**-3, orifice_length = 7*10**-3, manifold_length = 40*10**-3,
                sheath_l=25*10**-3)
    thermal = Thermal(engine.nozzle, k_x=237, k_y=237, material_Cp=1507.248, melting_point=933.5, material_density=2710)
                  
        