import numpy as np
from numpy import *
from constants import *
from engine import Engine

class Injector:
    
    def __init__(
        self,
        engine:Engine,
    ):
        self.engine = engine
        self.calculate()
        
        
    def calculate(self):
        # Determine pressure drop
        self.pressure_drop = 0.3 * self.engine.chamber_pressure
        print(self.pressure_drop)
        
        # Determine injection area
        self.injection_area = self.engine.ox_flow / (self.engine.discharge_coeff * sqrt(2 * self.pressure_drop * self.engine.ox_density))
        print(self.injection_area)
        
        # Determine injection velocity
        self.injection_velocity = sqrt(2 * self.pressure_drop / self.engine.ox_density)
        print(self.injection_velocity)
        
        # Determine manifold pressure
        self.pressure_manifold = self.engine.chamber_pressure + self.pressure_drop
        print(self.pressure_manifold)
        
