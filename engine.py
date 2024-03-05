from constants import *
from injector import Injector
from nozzle import Nozzle
from grain import Grain
import numpy as np
import matplotlib.pyplot as plt

class Engine:
    def __init__(
        self,
        Pc:float,
        Tc:float,
        thrust:float,
        M:float,
        mix_ratio:float,
        y:float,
        ox_den:float,
        cd:float,
        viscosity:float,
        orifice_diameter:float,
        tank_pressure:float,
        crit_pressure_drop:float,
        a:float,
        n:float,
        ox_flux:float,
        fuel_den:float,
        Pe:float = P_sl,
        Pa:float = P_sl
    ):
        print("")
        self.nozzle = Nozzle(Pc, Tc, thrust, M, mix_ratio, y)
        self.injector = Injector(self.nozzle, orifice_diam=orifice_diameter, tank_pressure=tank_pressure, crit_pressure_drop=crit_pressure_drop, Cd=cd, ox_den=ox_den, viscosity=viscosity)
        actual_thrust = self.nozzle.isp_m_s*(self.injector.ox_flow*(self.injector.mixture_ratio + 1)/self.injector.mixture_ratio)
        self.nozzle = Nozzle(Pc, Tc, actual_thrust, M, mix_ratio, y)
        self.grain = Grain(self.injector, a, n, fuel_den, ox_flux)
        
        self.describe()
        self.nozzle.plot()

        
    def describe(self):
        self.nozzle.describe()
        self.injector.describe()
        self.grain.describe()
        
    
        

engine = Engine(Pc=3*10**6, Tc=3391.91, thrust=300, M=0.026041, mix_ratio=5.3, y=1.2593, ox_den=786.6, cd=0.44, viscosity=3.237*10**-4,
                tank_pressure=50.525*10**5, crit_pressure_drop=17.4*10**5, orifice_diameter=1*10**-3, a=0.417, n=0.347, fuel_den=902, ox_flux=200)