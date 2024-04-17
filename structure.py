from nozzle import Nozzle
from grain import Grain
from injector import Injector
import numpy as np
import cadquery as cq

class Structure:
    def __init__(
        self,
        nozzle:Nozzle,
        grain:Grain,
        injector:Injector,
        bolt_OD:float,
        bolt_ID:float,
        bolt_distance:float,
        yield_stress:float,
        FoS:float,
        min_bolts:int = 4
    ):
        self.nozzle = nozzle
        self.grain=grain
        self.injector = injector
        self.bolt_OD = bolt_OD
        self.bolt_ID = bolt_ID
        self.bolt_distance = bolt_distance
        self.yield_stress = yield_stress
        self.FoS = FoS
        self.min_bolts = min_bolts
        
        self.calculate()
        
    def calculate(self):
        self.bolt_area = np.pi*self.bolt_ID**2 / 4
        self.tensile_force = np.pi * self.nozzle.radius_postcomb**2 * self.nozzle.chamber_pressure - self.nozzle.thrust
        self.bolts = int(np.ceil(self.FoS*(self.tensile_force/self.yield_stress)/self.bolt_area))
        if self.bolts < self.min_bolts:
            self.bolts = self.min_bolts
        self.FoS
        
    
    def model(self):
        # self.bolt_length = (
        #     self.grain.grain_length - self.nozzle.sheath_length*2 + self.nozzle.length + self.injector.precomb_length +self.injector.orifice_length
        #     + self.injector.manifold_length + 0.05
        # )*1000
        self.plate_width = (self.bolt_distance + 4*self.bolt_OD)*1000
        # self.bolt_geometry = (
        #     cq.Workplane("YZ").cylinder(self.bolt_length, 1000*self.bolt_OD/2)
        # )
        self.geometry = (
            cq.Workplane("XY").box(self.plate_width, self.plate_width, self.nozzle.plate_thickness*1000)
            .faces(">Z")
            .workplane()
            .hole((2*self.nozzle.clearance + self.nozzle.nozzle_OD - 2*self.nozzle.lip_thickness) * 1000)
            .faces(">Z")
            .polygon(self.bolts, self.bolt_distance*1000, forConstruction=True)
            .vertices()
            .hole((self.bolt_OD + 2*self.nozzle.clearance)*1000)
        )

    
    def describe(self):
        print("--------------STRUCTURE--------------")
        print(f"Tensile Force (N): {round(self.tensile_force)}")
        print("")