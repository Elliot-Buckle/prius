import numpy as np
from constants import *
from nozzle import Nozzle
import cadquery as cq
from ocp_vscode import *

class Injector:
    
    def __init__(
        self,
        nozzle:Nozzle,
        orifice_diam:float,
        tank_pressure:float,
        oxidizer:str,
        precomb_LD:float = 0.5,
        clearance:float = 0.1*10**-3,
        **kwargs
    ):
        self.nozzle = nozzle
        self.mixture_ratio = nozzle.mixture_ratio
        self.orifice_diameter = orifice_diam
        self.tank_pressure = tank_pressure
        self.critical_pressure_drop = 17.4*10**5
        self.discharge_coefficient = 0.77
        self.liquid_density = 786.5522078107688
        self.vapour_density = 158.08730384988698
        self.precomb_LD = precomb_LD
        self.clearance = clearance
        self.oxidizer = oxidizer
        self.calculate()
        
        # Checking if design dimensions were specified
        dimensions = True
        if "grain" in kwargs:
            self.grain = kwargs["grain"]
        else: dimensions = False
        
        if "lip_thickness" in kwargs:
            self.lip_thickness = kwargs["lip_thickness"]
        else: dimensions = False
        
        if "plate_thickness" in kwargs:
            self.plate_thickness = kwargs["plate_thickness"]
        else: dimensions = False
        
        if "sheath_length" in kwargs:
            self.sheath_length = kwargs["sheath_length"]
        else:
            dimensions = False
        if "orifice_length" in kwargs:
            self.orifice_length = kwargs["orifice_length"]
        else:
            dimensions = False
            
        if "injector_OD" in kwargs:
            self.injector_OD = kwargs["injector_OD"]
            self.injector_OR = self.injector_OD/2
        else:
            dimensions = False
            
        if "manifold_length" in kwargs:
            self.manifold_length = kwargs["manifold_length"]
        else:
            dimensions = False
            
        if dimensions is True:
            self.model()
        
    def describe(self):
        print("----------------INJECTOR----------------")
        print(f"Oxidizer mass flow rate (g/s): {round(self.ox_flow*10**3, 2)}")
        print(f"Pressure Drop (bar): {round(self.pressure_drop*10**-5, 2)}")
        print(f"Manifold Pressure (bar): {round(self.pressure_manifold*10**-5, 2)}")
        print(f"Number of orifices: {int(self.number_orifices)}")
        print(f"Injection velocity (m/s): {round(self.injection_velocity, 2)}")
        if self.oxidizer == "N2O":
            print(f"Condition: {self.condition}")
            print(f"Reynolds number: {round(self.reynolds_number, 2)}")
        print("")
        
    def model(self):
        self.precomb_diam = self.grain.outer_diameter - 2*self.lip_thickness
        self.precomb_radius = self.precomb_diam/2
        self.precomb_length = self.precomb_LD*self.precomb_diam
        height = self.precomb_length + self.sheath_length + self.orifice_length + self.manifold_length
        
        self.geometry = (
            cq.Workplane("YZ").cylinder(height*1000, self.injector_OR*1000)
            .faces(">X")
            .hole(self.grain.outer_diameter*1000 + 2*self.clearance*1000, self.sheath_length*1000)
            .faces(">X")
            .hole(self.precomb_diam*1000, self.precomb_length*1000 + self.sheath_length*1000)
            .faces("<X")
            .workplane()
            .hole(self.precomb_diam*1000, self.manifold_length*1000)
            .faces(">X")
            .workplane()
            .polygon(self.number_orifices, 1000*(self.precomb_diam/2),forConstruction=True,circumscribed=False)
            .vertices()
            .hole(self.orifice_diameter*1000)
            .faces("<X")
            .sketch()
            .circle(self.injector_OR*1000)
            .circle((self.injector_OR*1000 - self.lip_thickness*1000), mode="s")
            .finalize()
            .cutBlind(self.plate_thickness*1000)
        )
        cq.exporters.export(self.geometry, "injector.step")
    
    def calculate(self):
        #Determine oxidizer flow rate
        self.ox_flow = self.mixture_ratio*self.nozzle.thrust/(self.nozzle.isp_m_s*(self.mixture_ratio + 1))
        
        # Determine pressure drop
        self.pressure_drop = 0.3 * self.nozzle.chamber_pressure
        
        # Determine manifold pressure
        self.pressure_manifold = self.nozzle.chamber_pressure + self.pressure_drop
        
        # Determine orifice area
        self.orifice_area = 1/4*np.pi*self.orifice_diameter**2
        
        # Liquid Nitrous specific calcs
        if self.oxidizer=="N2O":
            if self.nozzle.chamber_pressure < self.tank_pressure - self.critical_pressure_drop:
                # Determine choked injection area
                self.injection_area = self.ox_flow / (self.discharge_coefficient * np.sqrt(2 * self.critical_pressure_drop * self.liquid_density))
                # Determine choked flow speed
                self.injection_velocity = np.sqrt(2 * self.critical_pressure_drop / self.liquid_density)
                self.condition = "choked"
                
            elif self.pressure_drop + self.nozzle.chamber_pressure > self.tank_pressure:
                # Caps manifold pressure at tank pressure
                self.pressure_drop = self.tank_pressure - self.nozzle.chamber_pressure
                # Determines injection area
                self.injection_area = self.ox_flow / (self.discharge_coefficient * np.sqrt(2 * self.pressure_drop * self.liquid_density))
                # Determine choked flow speed
                self.injection_velocity = np.sqrt(2 * self.pressure_drop / self.liquid_density)
                self.condition = "high chamber pressure"
            else:
                # Otherwise determine injection area
                self.injection_area = self.ox_flow / (self.discharge_coefficient * np.sqrt(2 * self.pressure_drop * self.liquid_density))
                # Determine flow speed
                self.injection_velocity = np.sqrt(2 * self.pressure_drop / self.liquid_density)
                self.condition = "single phase incompressible"
                
            
            # Determine number of orifices (rounded up)
            self.number_orifices = round(self.injection_area/self.orifice_area)
            
            # Determine actual oxidiser flow rate
            self.ox_flow = self.discharge_coefficient*self.injection_velocity*self.liquid_density*self.number_orifices*self.orifice_area

            #determine reynolds number
            self.reynolds_number = self.liquid_density*self.injection_velocity*self.orifice_diameter/(3.237*10**-4)
        # Gaseous oxygen specific calcs
        elif self.oxidizer=="GOX":
            self.number_orifices = round(426.1*self.ox_flow/(self.orifice_area*self.nozzle.chamber_pressure))
            self.injection_velocity = 165.8450071
        else:
            print("unrecognised oxidizer")
        
