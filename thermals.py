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
        
        nozzle.model(export=False, resolution=128)
        
        # Determine masses of cells in array
        
        # Create 2D array for cell masses
        grid_shape = np.shape(nozzle.grid_x)
        cell_array_shape = (grid_shape[0] - 1, grid_shape[1] - 1)
        print(grid_shape)
        cell_heat_capacities = np.full(cell_array_shape, NaN)
        
        # Create temperature array
        cell_temps = np.full(cell_array_shape, NaN)
        
        # Iterate over columns and rows
        for column in range(cell_array_shape[0]):
            for row in range(cell_array_shape[1]):
                # cell_volume = (nozzle.grid_x[cell - 1] + nozzle.grid_x[cell]) * ((nozzle.grid_y[cell] - nozzle.grid_y[cell - 1]) * pi)
                cell_area = abs(0.5*(nozzle.grid_x[row, column] - nozzle.grid_x[row, column + 1])*(nozzle.grid_y[row + 1, column] - nozzle.grid_y[row, column] + nozzle.grid_y[row, column + 1] - nozzle.grid_y[row, column]))
                cell_volume = abs(cell_area*np.pi*(nozzle.grid_y[row + 1, column] + nozzle.grid_y[row, column] + nozzle.grid_y[row, column + 1] + nozzle.grid_y[row, column])/2)
                #cell_volume = 0.5*(nozzle.grid_x_points[column + 1] - nozzle.grid_x_points[column])*((nozzle.grid_y_points[row + 1] - nozzle.grid_y_points[row + 1]) + (nozzle.grid_y_points[row] - nozzle.grid_y_points[row]))*2*pi*nozzle.grid_y_points[row]
                cell_mass = abs(cell_volume * self.material_rho)
                cell_heat_capacities[row, column] = cell_mass * self.material_Cp
        
        # Assigning temperature to cells
        cell_temps[np.nonzero(np.nan_to_num(cell_heat_capacities))] = 293.15
        
        # Calculating intercell distance (temporary approximation)
        intercell_dist = abs(nozzle.x_values[1] - nozzle.x_values[0])
        
        plt.imshow(cell_temps)
        plt.show()
                
                
        #plt.imshow(cell_heat_capacities)
        #plt.show()
        
        print("-----------------------")  
        #print(cell_heat_capacities)
        
        # Determine surface areas of cells in array
        # inner surface area
        cell_srf_areas_inner = np.full(cell_array_shape, NaN)
        # outer surface area
        cell_srf_areas_outer = np.full(cell_array_shape, NaN)
        # top surface area
        cell_srf_areas_aft = np.full(cell_array_shape, NaN)
        # bottom surface area
        cell_srf_areas_fore = np.full(cell_array_shape, NaN)
        
        for column in range(cell_array_shape[0]):
            for row in range(cell_array_shape[1]):
                # Inner surface area
                cell_srf_area = sqrt(
                    ((nozzle.grid_y[row, column] - nozzle.grid_y[row, column+1]) ** 2) + (nozzle.grid_x[row, column+1] - nozzle.grid_x[row, column]) ** 2 
                    ) * np.pi * (nozzle.grid_y[row, column] + nozzle.grid_y[row, column+1])
                cell_srf_areas_inner[row, column] = cell_srf_area
                
                # Outer surface area
                cell_srf_area = sqrt(
                    ((nozzle.grid_y[row + 1, column] - nozzle.grid_y[row + 1, column+1]) ** 2) + (nozzle.grid_x[row + 1, column+1] - nozzle.grid_x[row + 1, column]) ** 2 
                    ) * np.pi * (nozzle.grid_y[row + 1, column] + nozzle.grid_y[row + 1, column+1])
                cell_srf_areas_outer[row, column] = cell_srf_area
                
                # Aft surface area
                cell_srf_area = sqrt(
                    ((nozzle.grid_y[row, column + 1] - nozzle.grid_y[row + 1, column+1]) ** 2) + (nozzle.grid_x[row, column+1] - nozzle.grid_x[row + 1, column + 1]) ** 2 
                    ) * np.pi * (nozzle.grid_y[row + 1, column + 1] + nozzle.grid_y[row, column+1])
                cell_srf_areas_aft[row, column] = cell_srf_area
                
                # Fore surface area
                cell_srf_area = sqrt(
                    ((nozzle.grid_y[row, column] - nozzle.grid_y[row + 1, column]) ** 2) + (nozzle.grid_x[row, column] - nozzle.grid_x[row + 1, column]) ** 2 
                    ) * np.pi * (nozzle.grid_y[row + 1, column] + nozzle.grid_y[row, column])
                cell_srf_areas_fore[row, column] = cell_srf_area
                
        # plt.imshow(cell_srf_areas_fore)
        # #plt.imshow(cell_srf_areas_outer)
        # plt.show()
        
if __name__ == "__main__":
    #nozzle = Nozzle(Pc=3*10**6, Tc=3391.91, thrust=300, M=0.026041, mix_ratio=5.3, y=1.2593, nozzle_OD=75*10**-3)
    engine = Engine(fuel="HDPE", oxidizer="N2O", Pc=1.5*10**6, thrust=300, ox_den=786.6, cd=0.44,
                tank_pressure=50.525*10**5, crit_pressure_drop=17.4*10**5, orifice_diameter=1*10**-3, ox_flux=200,
                grain_OD=65*10**-3, cap_OD=75*10**-3, lip_t=10*10**-3, plate_t=10*10**-3, orifice_length = 10*10**-3, manifold_length = 20*10**-3,
                sheath_l=15*10**-3)
    thermal = Thermal(engine.nozzle, k_x=237, k_y=237, material_Cp=1507.248, melting_point=933.5, material_density=2710)
                  
        