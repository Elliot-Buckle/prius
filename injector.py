import numpy as np
from constants import *
from nozzle import Nozzle

class Injector:
    
    def __init__(
        self,
        nozzle:Nozzle,
        orifice_diam:float,
        tank_pressure:float,
        crit_pressure_drop:float,
        Cd:float,
        ox_den:float,
        viscosity:float
    ):
        self.nozzle = nozzle
        self.mixture_ratio = nozzle.mixture_ratio
        self.orifice_diameter = orifice_diam
        self.tank_pressure = tank_pressure
        self.critical_pressure_drop = crit_pressure_drop
        self.discharge_coefficient = Cd
        self.ox_density = ox_den
        self.dynamic_viscosity = viscosity
        
        self.calculate()
        
    def describe(self):
        print("----------------INJECTOR----------------")
        print(f"Oxidizer mass flow rate (g/s): {round(self.ox_flow*10**3, 2)}")
        print(f"Pressure Drop (bar): {round(self.pressure_drop*10**-5, 2)}")
        print(f"Manifold Pressure (bar): {round(self.pressure_manifold*10**-5, 2)}")
        print(f"Number of orifices: {int(self.number_orifices)}")
        print(f"Injection velocity (m/s): {round(self.injection_velocity, 2)}")
        print(f"Reynolds number: {round(self.reynolds_number, 2)}")
        print(f"Condition: {self.condition}")
        print("")
        
    def calculate(self):
        #Determine oxidizer flow rate
        self.ox_flow = self.mixture_ratio*self.nozzle.thrust/(self.nozzle.isp_m_s*(self.mixture_ratio + 1))
        
        # Determine pressure drop
        self.pressure_drop = 0.3 * self.nozzle.chamber_pressure
        if self.nozzle.chamber_pressure < self.tank_pressure - self.critical_pressure_drop:
            # Determine choked injection area
            self.injection_area = self.ox_flow / (self.discharge_coefficient * np.sqrt(2 * self.critical_pressure_drop * self.ox_density))
            # Determine choked flow speed
            self.injection_velocity = np.sqrt(2 * self.critical_pressure_drop / self.ox_density)
            self.condition = "choked"
            
        elif self.pressure_drop + self.nozzle.chamber_pressure > self.tank_pressure:
            # Caps manifold pressure at tank pressure
            self.pressure_drop = self.tank_pressure - self.nozzle.chamber_pressure
            # Determines injection area
            self.injection_area = self.ox_flow / (self.discharge_coefficient * np.sqrt(2 * self.pressure_drop * self.ox_density))
            # Determine choked flow speed
            self.injection_velocity = np.sqrt(2 * self.pressure_drop / self.ox_density)
            self.condition = "high chamber pressure"
        else:
            # Otherwise determine injection area
            self.injection_area = self.ox_flow / (self.discharge_coefficient * np.sqrt(2 * self.pressure_drop * self.ox_density))
            # Determine flow speed
            self.injection_velocity = np.sqrt(2 * self.pressure_drop / self.ox_density)
            self.condition = "single phase incompressible"
            
        
        # Determine orifice area
        self.orifice_area = 1/4*np.pi*self.orifice_diameter**2
        
        # Determine number of orifices (rounded up)
        self.number_orifices = np.ceil(self.injection_area/self.orifice_area)
        
        # Determine actual oxidiser flow rate
        self.ox_flow = self.discharge_coefficient*self.injection_velocity*self.ox_density*self.number_orifices*self.orifice_area
        
        # Determine manifold pressure
        self.pressure_manifold = self.nozzle.chamber_pressure + self.pressure_drop
        
        if self.dynamic_viscosity == 0:
            self.reynolds_number = 'nan'
        else:
            self.reynolds_number = self.ox_density*self.injection_velocity*self.orifice_diameter/self.dynamic_viscosity
        
