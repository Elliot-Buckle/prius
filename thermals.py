from constants import *
from nozzle import Nozzle
import numpy as np
from numpy import *
from ocp_vscode import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from rocketcea.cea_obj import CEA_Obj
from utilities import *
from materials import *

class Thermal:
    def __init__(
        self,
        nozzle:Nozzle,
        material:str,
        sim_file_name:str,
    ):
        self.conductivity_x = float(materials[material]["k_x"])
        self.conductivity_y = float(materials[material]["k_y"])
        self.material_Cp = float(materials[material]["C"])
        self.melting_point = float(materials[material]["Tmax"])
        self.material_rho = float(materials[material]["rho"])
        self.nozzle = nozzle
        self.material = material
        self.sim_file_name = sim_file_name
        self.calculate()
        
                
        # plt.imshow(self.cell_srf_areas_fore)
        # plt.imshow(self.cell_srf_areas_outer)
        # plt.show()
     
     
    # miscellaneous required calcs for heat sim           
    def calculate(self):
        self.nozzle.model()
        
                # Thickness in x direction
        cell_thickness = abs(self.nozzle.grid_x[-1, 0] - self.nozzle.grid_x[-1, 1])
        self.cell_thicknesses_x = cell_thickness
        
        # Thickness in y direction
        cell_thickness = abs(self.nozzle.grid_y[-1, 0] - self.nozzle.grid_y[-2, 0])
        self.cell_thicknesses_y = cell_thickness
        
        # Determine masses of cells in array
        
        # Create 2D array for cell masses
        grid_shape = np.shape(self.nozzle.grid_x)
        cell_array_shape = (grid_shape[0] - 1, grid_shape[1] - 1)
        self.cell_heat_capacities = np.full(cell_array_shape, NaN)
        
        # Create temperature array
        self.cell_temps = np.full(cell_array_shape, NaN)
        # Iterate over columns and rows
        # for column in range(cell_array_shape[0]):
        #     for row in range(cell_array_shape[1]):
        #         # cell_volume = (nozzle.grid_x[cell - 1] + nozzle.grid_x[cell]) * ((nozzle.grid_y[cell] - nozzle.grid_y[cell - 1]) * pi)
        #         cell_area = abs(0.5*(self.nozzle.grid_x[row, column] - self.nozzle.grid_x[row, column + 1])*(self.nozzle.grid_y[row + 1, column] - self.nozzle.grid_y[row, column] + self.nozzle.grid_y[row, column + 1] - self.nozzle.grid_y[row, column]))
        #         cell_volume = abs(cell_area*np.pi*(self.nozzle.grid_y[row, column] + self.nozzle.grid_y[row + 1, column] + self.nozzle.grid_y[row, column + 1] + self.nozzle.grid_y[row + 1, column + 1])/2)
        #         #cell_volume = 0.5*(nozzle.grid_x_points[column + 1] - nozzle.grid_x_points[column])*((nozzle.grid_y_points[row + 1] - nozzle.grid_y_points[row + 1]) + (nozzle.grid_y_points[row] - nozzle.grid_y_points[row]))*2*pi*nozzle.grid_y_points[row]
        #         cell_mass = abs(cell_volume * self.material_rho)
        #         self.cell_heat_capacities[row, column] = cell_mass * self.material_Cp
        
        self.cell_heat_capacities = self.material_Cp * self.material_rho * self.cell_thicknesses_x*self.cell_thicknesses_y* self.nozzle.grid_y[:-1, :-1]*2*np.pi
        # Assigning temperature to cells
        self.cell_temps[np.nonzero(np.nan_to_num(self.cell_heat_capacities))] = 293.15
                        
        #plt.imshow(self.cell_heat_capacities)
        #plt.show()
        
        print("-----------------------")  
        #print(self.cell_heat_capacities)
        
        # Determine surface areas and thicknesses of cells in array
        # inner surface area
        self.cell_srf_areas_inner = np.full(cell_array_shape, NaN)
        # outer surface area
        self.cell_srf_areas_outer = np.full(cell_array_shape, NaN)
        # top surface area
        self.cell_srf_areas_aft = np.full(cell_array_shape, NaN)
        # bottom surface area
        self.cell_srf_areas_fore = np.full(cell_array_shape, NaN)
        
        # # thickness in x-direction
        # self.cell_thicknesses_x = np.full(cell_array_shape, NaN)
        # # thickness in y-direction
        # self.cell_thicknesses_y = np.full(cell_array_shape, NaN)
        
        for column in range(cell_array_shape[0]):
            for row in range(cell_array_shape[1]):
                # Inner surface area
                cell_srf_area = sqrt(
                    ((self.nozzle.grid_y[row, column] - self.nozzle.grid_y[row, column+1]) ** 2) + (self.nozzle.grid_x[row, column+1] - self.nozzle.grid_x[row, column]) ** 2 
                    ) * np.pi * (self.nozzle.grid_y[row, column] + self.nozzle.grid_y[row, column+1])
                self.cell_srf_areas_inner[row, column] = cell_srf_area
                
                # Outer surface area
                cell_srf_area = sqrt(
                    ((self.nozzle.grid_y[row + 1, column] - self.nozzle.grid_y[row + 1, column+1]) ** 2) + (self.nozzle.grid_x[row + 1, column+1] - self.nozzle.grid_x[row + 1, column]) ** 2 
                    ) * np.pi * (self.nozzle.grid_y[row + 1, column] + self.nozzle.grid_y[row + 1, column+1])
                self.cell_srf_areas_outer[row, column] = cell_srf_area
                
                # Aft surface area
                cell_srf_area = sqrt(
                    ((self.nozzle.grid_y[row, column + 1] - self.nozzle.grid_y[row + 1, column+1]) ** 2) + (self.nozzle.grid_x[row, column+1] - self.nozzle.grid_x[row + 1, column + 1]) ** 2 
                    ) * np.pi * (self.nozzle.grid_y[row + 1, column + 1] + self.nozzle.grid_y[row, column+1])
                self.cell_srf_areas_aft[row, column] = cell_srf_area
                
                # Fore surface area
                cell_srf_area = sqrt(
                    ((self.nozzle.grid_y[row, column] - self.nozzle.grid_y[row + 1, column]) ** 2) + (self.nozzle.grid_x[row, column] - self.nozzle.grid_x[row + 1, column]) ** 2 
                    ) * np.pi * (self.nozzle.grid_y[row + 1, column] + self.nozzle.grid_y[row, column])
                self.cell_srf_areas_fore[row, column] = cell_srf_area
                
        # self.cell_heat_capacities = self.material_Cp * self.material_rho * self.cell_thicknesses_x*self.cell_thicknesses_y* self.nozzle.grid_y[:-1, :-1]
        # Calculate gas properties
        self.mol_wt_gamma = self.nozzle.prop.get_Chamber_MolWt_gamma(self.nozzle.Pc_psi, self.nozzle.mixture_ratio)
        self.gamma = self.mol_wt_gamma[1]
        self.mol_wt = self.mol_wt_gamma[0]/1000 # kg/mol
        # Index of throat values
        self.throat_index = np.argmin(self.nozzle.gas_side_y_values)
        # Subsonic machs
        subsonic_mach = mach_from_area(np.pi*self.nozzle.gas_side_y_values[:self.throat_index]**2, self.nozzle.throat_area, self.gamma, is_subsonic=True)
        # Supersonic machs
        supersonic_mach = mach_from_area(np.pi*self.nozzle.gas_side_y_values[self.throat_index:]**2, self.nozzle.throat_area, self.gamma, is_subsonic=False)
        # array of machs for cells
        self.mach_numbers = np.concatenate((subsonic_mach, supersonic_mach))
        # calculate temperatures
        self.chamber_temp = self.nozzle.prop.get_Tcomb(self.nozzle.Pc_psi, self.nozzle.mixture_ratio)*5/9 #Kelvin
        self.gas_temperatures = temp_from_mach(self.mach_numbers, self.chamber_temp, self.gamma)
        # calculate mass fluxes
        self.mass_fluxes = mass_flux(self.mach_numbers, self.gamma, self.nozzle.chamber_pressure, self.chamber_temp, self.mol_wt)
        # get variables (who decided that cal/g K was a good unit for Cp smh)
        
        # returns tuple of Cp, viscosity, conductivity, prandtl number
        chamber_transport_properties = self.nozzle.prop.get_Chamber_Transport(self.nozzle.Pc_psi, self.nozzle.mixture_ratio, frozen=1)
        throat_transport_properties = self.nozzle.prop.get_Throat_Transport(self.nozzle.Pc_psi, self.nozzle.mixture_ratio, frozen=1)
        exit_transport_properties = self.nozzle.prop.get_Exit_Transport(self.nozzle.Pc_psi, self.nozzle.mixture_ratio, frozen=1)

        # Interpolating transport properties
        # array of chamber, throat and exit temperatures
        CTE_temps = np.array([self.gas_temperatures[0], self.gas_temperatures[self.throat_index], self.gas_temperatures[-1]])
        # array of chamber, throat, exit prandtl numbers
        CTE_Pr = np.array([chamber_transport_properties[3], throat_transport_properties[3], exit_transport_properties[3]])
        # Interpolated prandtl number
        self.prandtl_numbers = np.flip(np.interp(np.flip(self.gas_temperatures), np.flip(CTE_temps) , np.flip(CTE_Pr)))
        
        # array of chamber, throat, exit conductivities
        CTE_k = np.array([chamber_transport_properties[2], throat_transport_properties[2], exit_transport_properties[2]])
        # Interpolated conductivites
        self.gas_conductivities = np.flip(np.interp(np.flip(self.gas_temperatures), np.flip(CTE_temps) , np.flip(CTE_k)))
        
        # array of chamber, throat, exit viscosities (Pa s)
        CTE_mu = np.array([chamber_transport_properties[1], throat_transport_properties[1], exit_transport_properties[1]])*0.0001
        # Interpolated viscosity
        self.gas_viscosities = np.flip(np.interp(np.flip(self.gas_temperatures), np.flip(CTE_temps) , np.flip(CTE_mu)))
        
        #transport_properties = ispObj.get_Chamber_Transport(Pc=Nozzle.Pc_psi, MR=Nozzle.mixture_ratio, frozen=1)
     
    # Determine heat flux on chamber wall
    def wall_heat_fluxes(self, wall_temp):
        wall_heat_flux = 0.023 * self.mass_fluxes[:-1] ** 0.8 * self.prandtl_numbers[:-1] ** 0.4 * self.gas_conductivities[:-1] / (self.nozzle.gas_side_y_values[:-1]**0.2 * self.gas_viscosities[:-1]**0.8)
        #print(((self.nozzle.gas_side_y_values[:-1] ** 0.2) * (self.gas_conductivities[:-1] ** 0.8))) * (self.chamber_temp - wall_temp
        return wall_heat_flux
    
    # Determine heat flux on adjacent cells
    def cell_heat_fluxes(self, material_conductivity, cell_thickness, cell_temp, adjacent_cell_temp):
        cell_heat_flux = (material_conductivity / cell_thickness) * (adjacent_cell_temp - cell_temp)
        return cell_heat_flux
    
    # Main simulation
    def simulation(self):
        timestep = 0.0001
        self.nozzle.model()
        i=0
        time = 0
        # Create plots for animation
        fig, ax = plt.subplots()
        ims = []
        
        # Main loop
        while np.nanmax(self.cell_temps) < self.melting_point:
            # Create array for adjacent cell temperature
            # Apply heat flux to cells along wall
            self.wall_delta_heat = timestep * self.cell_srf_areas_outer[self.nozzle.cell_wall_rows, self.nozzle.cell_wall_columns] * self.wall_heat_fluxes(self.cell_temps[self.nozzle.cell_wall_rows, self.nozzle.cell_wall_columns])
            delta_T = np.nan_to_num(self.wall_delta_heat / self.cell_heat_capacities[self.nozzle.cell_wall_rows, self.nozzle.cell_wall_columns])
            self.cell_temps[self.nozzle.cell_wall_rows, self.nozzle.cell_wall_columns] += delta_T
            # Apply heat flux to cells below
            roll_down = np.roll(self.cell_temps, 1, 0)
            roll_down[0, :] = NaN
            roll_up = np.roll(self.cell_temps, -1, 0)
            roll_up[-1, :] = self.cell_temps[-1, :]
            roll_right = np.roll(self.cell_temps, 1, 1)
            roll_right[:, 0] = self.cell_temps[:, 0]
            roll_left = np.roll(self.cell_temps, -1, 1)
            roll_left[:, -1] = self.cell_temps[:, -1]
            
            self.delta_cell_heat = timestep * (
                self.cell_srf_areas_outer * self.cell_heat_fluxes(self.conductivity_y, self.cell_thicknesses_y, self.cell_temps, roll_down)
            +self.cell_srf_areas_inner * self.cell_heat_fluxes(self.conductivity_y, self.cell_thicknesses_y, self.cell_temps, roll_up)
            + self.cell_srf_areas_aft * self.cell_heat_fluxes(self.conductivity_x, self.cell_thicknesses_x, self.cell_temps, roll_right)
            + self.cell_srf_areas_fore * self.cell_heat_fluxes(self.conductivity_x, self.cell_thicknesses_x, self.cell_temps, roll_left)
            )
            delta_T = np.nan_to_num(self.delta_cell_heat / self.cell_heat_capacities)
            self.cell_temps += delta_T
            
            # Create frames for animation and add title with timestamp
            if i % int(round(1/(10*timestep))) == 0:
                im = ax.imshow(self.cell_temps, cmap="inferno", animated=True)
                ann = ax.annotate(f"{self.material} Nozzle, t = {round(time)}s", (0.5, 1.03), xycoords="axes fraction", ha="center")
                ims.append([im, ann])
            
            time += timestep
            i+=1
        self.max_burn_time = time
        
        print("--------------------THERMAL--------------------")
        print(f"Time before safe temp limit exceeded: {round(self.max_burn_time, 1)}s")
        print(f"Maximum Impulse (Ns): {round(self.nozzle.thrust*self.max_burn_time)}")
        print("")
        # Add colour bar to animation and save file
        fig.colorbar(mappable=im, ax=ax, label="Temperature (K)")
        sim = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=1000)
        sim.save(self.sim_file_name)
        

        
            
            
            
            
            
        
if __name__ == "__main__":
    #nozzle = Nozzle(Pc=3*10**6, Tc=3391.91, thrust=300, M=0.026041, mix_ratio=5.3, y=1.2593, nozzle_OD=75*10**-3)
    pass