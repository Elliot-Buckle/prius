from constants import *
from injector import Injector
from nozzle import Nozzle
from grain import Grain
from structure import Structure
import numpy as np
import matplotlib.pyplot as plt
import cadquery as cq
from ocp_vscode import *
from rocketcea.cea_obj import CEA_Obj
from propellants import *
from thermals import Thermal
from materials import *

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
        material:str = "Aluminium",
        sim_file_name:str = "heat_sim",
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
        self.material = material
        self.sim_file_name = sim_file_name+".gif"
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
        
        # Generating Structure
        self.structure = Strucutre(self.nozzle, bolt_OD=5*10**-3, bolt_ID=4.134*10**-3, bolt_distance=100*10**-3, yield_stress=400*10**6, FoS=4, grain=self.grain, injector=self.injector)
        self.thermal = Thermal(self.nozzle, self.material, self.sim_file_name)
        self.model()
        #show(self.grain.geometry)
        #self.injector.model()
        #self.nozzle.plot()
        #show(self.nozzle.geometry)
        self.describe()
        #self.nozzle.plot()

    def describe(self):
        print("")
        print("")
        print("")
        print("")
        print("")
        print("----------------ENGINE----------------")
        print("Fuel: " + self.fuel)
        print("Oxidizer: " + self.oxidizer)
        print(f"O/F ratio: {round(self.nozzle.mixture_ratio, 2)}")
        print("")
        self.grain.describe()
        self.nozzle.describe()
        self.injector.describe()
        self.structure.describe()
        
    def model(self, export=False):
        grain_colour = colours[self.fuel]
        nozzle_colour = materials[self.material]["colour"]
        self.grain.model()
        self.injector.model()
        self.nozzle.model()
        self.structure.model()
        self.engine = cq.Assembly(name="Engine")
        # self.engine.add(self.nozzle.geometry, name="nozzle", loc=cq.Location((100,0,0), (1,0,0), 0))
        self.engine.add(self.nozzle.geometry, name="nozzle", color=nozzle_colour)
        self.engine.add(self.grain.geometry, name="grain", color=grain_colour)
        self.engine.add(self.injector.geometry, name="injector", color="silver")
        self.engine.add(self.structure.geometry, name="plate1", color="gray")
        self.engine.add(self.structure.geometry, name="plate2", color="gray")
        # for i in range(self.structure.bolts):
        #     self.engine.add(self.structure.bolt_geometry, name="bolt"+str(i+1), color="gray")
        # aligning bolt holes
        self.engine.constrain("plate2@faces@>X", "plate1@faces@>X", "Axis")
        # aligning nozzle plate
        self.engine.constrain("plate1@faces@>Z", "nozzle@faces@>X", "Axis")
        self.engine.constrain("plate1@faces@<Z", "nozzle@faces@>X", "Point")
        # aligning injector plate
        self.engine.constrain("plate2@faces@>Z", "injector@faces@>X", "Axis")
        self.engine.constrain("plate2@faces@>Z", "injector@faces@<X", "Point")
        # aligning faces of nozzle and grain
        self.engine.constrain("nozzle@faces@>Z", "grain@faces@<X", "Axis")
        # putting two together
        self.engine.constrain("nozzle@faces@>X[1]", "grain@faces@<X", "Point")
        # aligning faces of injector and grain
        self.engine.constrain("injector@faces@<X[1]", "grain@faces@>X", "Plane")
        self.engine.solve()
        
        show(self.engine)
        
    
        
# Oxidizer can be either "N2O" or "GOX"
engine = Engine(fuel="HDPE", oxidizer="GOX", Pc=3*10**6, thrust=300, ox_den=786.6, cd=0.44,
                tank_pressure=50.525*10**5, crit_pressure_drop=17.4*10**5, orifice_diameter=1*10**-3, ox_flux=200,
                grain_OD=65*10**-3, cap_OD=75*10**-3, lip_t=10*10**-3, plate_t=10*10**-3, orifice_length = 10*10**-3, manifold_length = 20*10**-3,
                sheath_l=15*10**-3, material="Aluminium", sim_file_name="Ally")
engine.thermal.simulation()