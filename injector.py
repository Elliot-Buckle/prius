import numpy as np
from numpy import *
from constants import *
from engine import Engine

class Injector:
    
    def __init__(
        self,
        engine:Engine,
    ):
    
        self.calculate()
        
        
    def calculate(self):
        # Determine pressure drop
        self.pressure_drop = 0.3 * Engine.pressure_chamber
        print(self.pressure_drop)
        
        # Determine injection area
        self.injection_area = Engine.ox_flow / (Engine.discharge_coefficient * sqrt(2 * self.pressure_drop * Engine.ox_density))
        print(self.injection_area)
        
        # Determine injection velocity
        self.injection_velocity = sqrt(2 * self.pressure_drop / Engine.ox_density)
        print(self.injection_velocity)
        
        # Determine manifold pressure
        self.pressure_manifold = Engine.pressure_chamber + self.pressure_drop
        print(self.pressure_manifold)
        
