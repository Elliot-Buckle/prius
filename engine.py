from constants import *
from injector import Injector
from nozzle import Nozzle
from grain import Grain
import numpy as np
import matplotlib.pyplot as plt
import cadquery as cq
from ocp_vscode import *

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
        grain_OD:float,
        nozzle_OD:float,
        postcomb_LD:float = 1,
        lip_t:float = 5*10**-3,
        plate_t:float = 5*10**-3,
        Pe:float = P_sl,
        Pa:float = P_sl
    ):
        print("")
        self.nozzle = Nozzle(Pc, Tc, thrust, M, mix_ratio, y, nozzle_OD)
        self.injector = Injector(self.nozzle, orifice_diam=orifice_diameter, tank_pressure=tank_pressure, crit_pressure_drop=crit_pressure_drop, Cd=cd, ox_den=ox_den, viscosity=viscosity)
        actual_thrust = self.nozzle.isp_m_s*(self.injector.ox_flow*(self.injector.mixture_ratio + 1)/self.injector.mixture_ratio)
        self.grain = Grain(self.injector, a, n, fuel_den, ox_flux, grain_OD)
        self.nozzle = Nozzle(Pc, Tc, actual_thrust, M, mix_ratio, y, nozzle_OD, plate_thickness=plate_t, lip_thickness=lip_t, grain=self.grain,sheath_length=25*10**-3)
        self.grain.model()
        #self.nozzle.plot()
        show(self.nozzle.geometry)
        self.describe()
        # self.nozzle.plot()

    def describe(self):
        self.nozzle.describe()
        self.injector.describe()
        self.grain.describe()
        
    
        

engine = Engine(Pc=3*10**6, Tc=3391.91, thrust=300, M=0.026041, mix_ratio=5.3, y=1.2593, ox_den=786.6, cd=0.44, viscosity=3.237*10**-4,
                tank_pressure=50.525*10**5, crit_pressure_drop=17.4*10**5, orifice_diameter=1*10**-3, a=0.417, n=0.347, fuel_den=902, ox_flux=200,
                grain_OD=65*10**-3, nozzle_OD=85*10**-3, lip_t=10*10**-3, plate_t=10*10**-3)