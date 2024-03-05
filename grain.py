from injector import Injector
import numpy as np
class Grain:
    def __init__(
        self,
        injector:Injector,
        a:float,
        n:float,
        fuel_den:float,
        ox_flux:float,
        geometry:str = "tube"
    ):
        self.injector = injector
        self.regression_constant = a
        self.regression_exponent = n
        self.ox_flux = ox_flux
        self.geometry = geometry.casefold()
        self.fuel_density = fuel_den
        
        self.calculate()
        
    def calculate(self):
        self.regression_rate = 10**-3*self.regression_constant*(self.ox_flux/10)**self.regression_exponent
        
        if self.geometry == "tube":
            self.port_area = self.injector.ox_flow/self.ox_flux
            self.port_radius = np.sqrt(self.port_area/np.pi)
            self.port_diameter = 2*self.port_radius
            self.port_circumference = 2*np.pi*self.port_radius
            self.fuel_flow_rate = self.injector.ox_flow/self.injector.mixture_ratio
            self.grain_length = self.fuel_flow_rate/(self.regression_rate*self.fuel_density*self.port_circumference)
            
    def describe(self):
        print("----------------GRAIN----------------")
        print(f"Port Diameter (mm): {round(self.port_diameter*10**3, 2)}")
        print(f"Grain Length (mm): {round(self.grain_length*10**3, 2)}")
        print(f"Fuel mass flowrate (g/s): {round(self.fuel_flow_rate*10**3, 2)}")
        print("")