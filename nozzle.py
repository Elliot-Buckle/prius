from constants import *
#from injector import Injector
from numpy import *
import numpy as np
from math import *
import matplotlib.pyplot as plt


class Nozzle:
    def __init__(
        self,
        Pc:float,
        Tc:float,
        thrust:float,
        M:float,
        mix_ratio:float,
        y:float,
        Pe:float = P_sl,
        Pa:float = P_sl,
        contour:str = "moc",
        N:int = 2000
    ):
        self.chamber_pressure = Pc
        self.exit_pressure_1D = Pe
        self.chamber_temp = Tc
        self.mixture_ratio = mix_ratio
        self.ratio_of_heats = y
        self.mol_mass = M
        self.thrust = thrust
        self.ambient_pressure = Pa
        self.contour = contour.casefold()
        self.xpoints = None
        self.ypoints = None
        self.divisions = N
        
        self.calculate()
        self.plot()
        
    def calculate(self):
        
        self.throat_pressure = self.chamber_pressure*((self.ratio_of_heats + 1)/2)**(self.ratio_of_heats/(1 - self.ratio_of_heats))
        
        self.area_ratio_1D = 1/(((self.ratio_of_heats + 1)/2)**(1/(self.ratio_of_heats - 1))*(self.exit_pressure_1D/self.chamber_pressure)**(1/self.ratio_of_heats)
                              *np.sqrt((self.ratio_of_heats + 1)/(self.ratio_of_heats - 1)*(1 - (self.exit_pressure_1D/self.chamber_pressure)**((self.ratio_of_heats - 1)/self.ratio_of_heats)))
                              )
        
        self.exit_vel_1D = np.sqrt(
            2*self.ratio_of_heats/
            (self.mol_mass*(self.ratio_of_heats - 1))*R*self.chamber_temp*(1 - (self.exit_pressure_1D/self.chamber_pressure)**((self.ratio_of_heats - 1)/self.ratio_of_heats))
            )
        
        self.throat_mass_flux = self.chamber_pressure*self.ratio_of_heats*np.sqrt(
            (2/(self.ratio_of_heats + 1))**((self.ratio_of_heats + 1)/(self.ratio_of_heats - 1))
        )/np.sqrt(self.ratio_of_heats*R*self.chamber_temp)
        
        self.isp_m_s = (self.exit_vel_1D*self.throat_mass_flux + self.area_ratio_1D*(self.exit_pressure_1D - self.ambient_pressure))/self.throat_mass_flux
        
        self.isp_s = self.isp_m_s/g
        
        self.mass_flow = self.thrust/self.isp_m_s
        
        self.throat_area = self.mass_flow/self.throat_mass_flux
        
        self.throat_radius = np.sqrt(self.throat_area/np.pi)
        self.throat_diameter = 2*self.throat_radius
        self.exit_area_1D = self.throat_area * self.area_ratio_1D
        self.exit_radius_1D = np.sqrt(self.exit_area_1D/np.pi)
        self.exit_diameter_1D = 2*self.exit_radius_1D
        
    def describe(self):
        print("----------------NOZZLE----------------")
        print(f"Thrust (N): {round(self.thrust,2)}")
        print(f"Isp (s): {round(self.isp_s, 2)}")
        print(f"Isp (Ns/kg): {round(self.isp_m_s, 0)}")
        print(f"Throat Diameter (mm): {round(self.throat_diameter*10**3, 2)}")
        print(f"Exit Diameter (mm): {round(self.exit_diameter_1D*10**3, 2)}")
        print(f"Expansion ratio: {round(self.area_ratio_1D, 2)}")
        print("")
        
    def generate_contour(self):
        # temp rao nozzle code
        
        # define some required variables
        length_fraction = 0.8
        parabola_angle_initial = 22
        parabola_angle_final = 13
        
        # determine length of nozzle
        self.length = length_fraction * (
            ((sqrt(self.expansion_ratio) - 1) * self.radius_throat)
            / (tan(radians(self.angle_divergent)))
        )
        
        # converging throat section
        t1  = linspace(radians(-135), radians(-90), self.divisions, False)
        
        x_conv = 1.5 * self.throat_radius * cos(t1)
        y_conv = (1.5 * self.throat_radius * sin(t1)) + (2.5 * self.throat_radius)
        
        # diverging throat section
        t2 = linspace(radians(-90), radians(parabola_angle_initial - 90), self.divisions, False)
        
        x_div = 0.382 * self.throat_radius * cos(t2)
        y_div = (0.382 * self.throat_radius * sin(t2)) + (1.382 * self.throat_radius)
        
        # start of parabola
        Nx = (
            0.382
            * self.throat_radius
            * cos(radians(parabola_angle_initial) - radians(90))
        )
        
        Ny = (
            0.382
            * self.throat_radius
            * sin(radians(parabola_angle_initial) - radians(90))
        ) + (1.382 * self.throat_radius)   
        
        # intersection of line segments
        Qx = (
            (self.exit_radius_1D - (self.length * tan(radians(parabola_angle_final))))
            - (Ny - (Nx * tan(radians(parabola_angle_initial))))
        ) / (
            tan(radians(parabola_angle_initial))
            - tan(radians(parabola_angle_final))
        )
        
        Qy = (
            (
                (
                    tan(radians(parabola_angle_initial))
                    * (
                        self.exit_radius_1D
                        - (self.length * tan(radians(parabola_angle_final)))
                    )
                )
            )
            - (
                tan(radians(parabola_angle_final))
                * (Ny - (Nx * tan(radians(parabola_angle_initial))))
            )
        ) / (
            tan(radians(parabola_angle_initial))
            - tan(radians(parabola_angle_final))
        )
        
        # parabolic expansion section
        t3 = linspace(0, 1, self.divisions, False)
        
        x_para = (
            (Nx * ((1 - t3) ** 2))
            + (Qx * (2 * t3) * (1 - t3))
            + ((t3**2) * self.length)
        )
        y_para = (
            (Ny * ((1 - t3) ** 2))
            + (Qy * (2 * t3) * (1 - t3))
            + ((t3**2) * self.radius_exit)
        )
        
        # add contour to x/y points
        self.xpoints = concatenate((x_conv, x_div, x_para))
        self.ypoints = concatenate((y_conv, y_div, y_para))
        
        
        
#nozzle = Nozzle(Pc=2*10**6, Tc=3000, thrust=300, M=0.026, mix_ratio=8, y=1.2, ox_den=700)