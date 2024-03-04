import numpy as np
from constants import *
from nozzle import Nozzle

class Injector:
    
    def __init__(
        self,
        nozzle:Nozzle,
    ):
        self.nozzle = nozzle
        self.calculate()
        
        
    def calculate(self):
        # Determine pressure drop
        self.pressure_drop = 0.3 * self.nozzle.chamber_pressure
        print(self.pressure_drop)
        
        # Determine injection area
        self.injection_area = self.nozzle.ox_flow / (self.nozzle.discharge_coeff * np.sqrt(2 * self.pressure_drop * self.nozzle.ox_density))
        print(self.injection_area)
        
        # Determine injection velocity
        self.injection_velocity = np.sqrt(2 * self.pressure_drop / self.nozzle.ox_density)
        print(self.injection_velocity)
        
        # Determine manifold pressure
        self.pressure_manifold = self.nozzle.chamber_pressure + self.pressure_drop
        print(self.pressure_manifold)
        
