from constants import *
from injector import Injector
from nozzle import Nozzle
from grain import Grain
import numpy as np
import matplotlib.pyplot as plt
import cadquery as cq
from ocp_vscode import *
from rocketcea.cea_obj import CEA_Obj

class Engine:
    def __init__(
        self,
        fuel:str,
        oxidizer:str,
        Pc:float,
        thrust:float,
        orifice_diameter:float,
        tank_pressure:float,
        ox_flux:float,
        grain_OD:float,
        cap_OD:float,
        postcomb_LD:float = 1,
        lip_t:float = 5*10**-3,
        plate_t:float = 5*10**-3,
        Pe:float = P_sl,
        Pa:float = P_sl,
        **kwargs
    ):
        if "sheath_l" in kwargs:
            sheath_l = kwargs["sheath_l"]
        else:
            sheath_l = 0
            
        self.fuel = fuel
        self.oxidizer = oxidizer
        print("")
        # Initial call of nozzle to calculate approximate flow rates
        self.nozzle = Nozzle(fuel=fuel, oxidizer=oxidizer, Pc=Pc, thrust=thrust, nozzle_OD=cap_OD, Pe=Pe, Pa=Pa)
        
        # Initial call of injector to calculate orifice no. and actual flow rate
        self.injector = Injector(nozzle=self.nozzle, orifice_diam=orifice_diameter, tank_pressure=tank_pressure, oxidizer=oxidizer)
        actual_thrust = self.nozzle.isp_m_s*(self.injector.ox_flow*(self.injector.mixture_ratio + 1)/self.injector.mixture_ratio)
        
        # Call of grain
        self.grain = Grain(self.injector, fuel, oxidizer, ox_flux, grain_OD)
        
        # Second call of injector to do modelling based on grain results
        self.injector = Injector(self.nozzle, orifice_diam=orifice_diameter, tank_pressure=tank_pressure, grain=self.grain, lip_thickness=lip_t, plate_thickness=plate_t,
                                 sheath_length=kwargs["sheath_l"], manifold_length=kwargs["manifold_length"], orifice_length=kwargs["orifice_length"], injector_OD=cap_OD, oxidizer=oxidizer)
        
        # Second call of nozzle to do modelling based on grain results
        self.nozzle = Nozzle(fuel, oxidizer, Pc, actual_thrust, cap_OD, plate_thickness=plate_t, lip_thickness=lip_t, grain=self.grain, sheath_length=sheath_l, Pe=Pe, Pa=Pa)
        
        self.model(export=True)
        #show(self.grain.geometry)
        #self.injector.model()
        #self.nozzle.plot()
        #show(self.nozzle.geometry)
        self.describe()
        #self.nozzle.plot()

    def describe(self):
        print("----------------ENGINE----------------")
        print("Fuel: " + self.fuel)
        print("Oxidizer: " + self.oxidizer)
        print(f"O/F ratio: {round(self.nozzle.mixture_ratio, 2)}")
        print("")
        self.nozzle.describe()
        self.injector.describe()
        self.grain.describe()
        
    def model(self, export=False):
        self.grain.model()
        self.injector.model()
        self.nozzle.model()
        self.engine = cq.Assembly(name="Engine")
        # self.engine.add(self.nozzle.geometry, name="nozzle", loc=cq.Location((100,0,0), (1,0,0), 0))
        self.engine.add(self.nozzle.geometry, name="nozzle", color="gray")
        self.engine.add(self.grain.geometry, name="grain", color="white")
        self.engine.add(self.injector.geometry, name="injector", color="gray")
        # aligning faces of nozzle and grain
        self.engine.constrain("nozzle@faces@>Z", "grain@faces@<X", "Axis")
        # putting two together
        self.engine.constrain("nozzle@faces@>X[1]", "grain@faces@<X", "Point")
        # aligning faces of injector and grain
        self.engine.constrain("injector@faces@<X[1]", "grain@faces@>X", "Plane")
        self.engine.solve()
        
        show(self.engine)
        
    
        
# Oxidizer can be either N2O or GOX
engine = Engine(fuel="HDPE", oxidizer="N2O", Pc=1.5*10**6, thrust=300, ox_den=786.6, cd=0.44,
                tank_pressure=50.525*10**5, crit_pressure_drop=17.4*10**5, orifice_diameter=1*10**-3, ox_flux=200,
                grain_OD=65*10**-3, cap_OD=75*10**-3, lip_t=10*10**-3, plate_t=10*10**-3, orifice_length = 10*10**-3, manifold_length = 20*10**-3,
                sheath_l=15*10**-3)