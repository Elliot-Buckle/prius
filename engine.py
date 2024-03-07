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
        cap_OD:float,
        postcomb_LD:float = 1,
        lip_t:float = 5*10**-3,
        plate_t:float = 5*10**-3,
        Pe:float = P_sl,
        Pa:float = P_sl,
        **kwargs
    ):
        print("")
        # Initial call of nozzle to calculate approximate flow rates
        self.nozzle = Nozzle(Pc, Tc, thrust, M, mix_ratio, y, cap_OD)
        
        # Initial call of injector to calculate orifice no. and actual flow rate
        self.injector = Injector(self.nozzle, orifice_diam=orifice_diameter, tank_pressure=tank_pressure, crit_pressure_drop=crit_pressure_drop, Cd=cd, ox_den=ox_den, viscosity=viscosity)
        actual_thrust = self.nozzle.isp_m_s*(self.injector.ox_flow*(self.injector.mixture_ratio + 1)/self.injector.mixture_ratio)
        
        # Call of grain
        self.grain = Grain(self.injector, a, n, fuel_den, ox_flux, grain_OD)
        
        # Second call of injector to do modelling based on grain results
        self.injector = Injector(self.nozzle, orifice_diam=orifice_diameter, tank_pressure=tank_pressure, crit_pressure_drop=crit_pressure_drop, Cd=cd, ox_den=ox_den, viscosity=viscosity, grain=self.grain, lip_thickness=lip_t, plate_thickness=plate_t,
                                 sheath_length=kwargs["sheath_l"], manifold_length=kwargs["manifold_length"], orifice_length=kwargs["orifice_length"], injector_OD=cap_OD)
        
        # Second call of nozzle to do modelling based on grain results
        self.nozzle = Nozzle(Pc, Tc, actual_thrust, M, mix_ratio, y, cap_OD, plate_thickness=plate_t, lip_thickness=lip_t, grain=self.grain, sheath_length=25*10**-3)
        
        self.model()
        show(self.nozzle.geometry)
        #self.injector.model()
        #self.nozzle.plot()
        #show(self.nozzle.geometry)
        self.describe()
        # self.nozzle.plot()

    def describe(self):
        self.nozzle.describe()
        self.injector.describe()
        self.grain.describe()
        
    def model(self):
        self.grain.model()
        self.injector.model()
        self.nozzle.model()
        self.engine = cq.Assembly()
        self.engine.add(self.nozzle.geometry, name="nozzle")
        self.engine.add(self.grain.geometry, name="grain")
        self.engine.add(self.injector.geometry, name="injector")
        #self.engine.constrain("grain@faces@>X", )
        
    
        

engine = Engine(Pc=3*10**6, Tc=3252.15, thrust=300, M=0.025423, mix_ratio=6.1, y=1.2501, ox_den=786.6, cd=0.44, viscosity=3.237*10**-4,
                tank_pressure=50.525*10**5, crit_pressure_drop=17.4*10**5, orifice_diameter=1*10**-3, a=0.417, n=0.347, fuel_den=902, ox_flux=200,
                grain_OD=65*10**-3, cap_OD=75*10**-3, lip_t=10*10**-3, plate_t=10*10**-3, orifice_length = 7*10**-3, manifold_length = 40*10**-3,
                sheath_l=25*10**-3)