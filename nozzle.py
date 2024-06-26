from constants import *
#from injector import Injector
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
import cadquery as cq
from ocp_vscode import *
from propellants import *
from rocketcea.cea_obj import CEA_Obj, add_new_fuel, add_new_oxidizer
#from grain import Grain
    

class Nozzle:
    def __init__(
        self,
        fuel:str,
        oxidizer:str,
        Pc:float,
        thrust:float,
        nozzle_OD:float,
        angle_converging:float = np.pi/3,
        postcomb_LD:float = 1,
        Pe:float = P_sl,
        Pa:float = P_sl,
        contour:str = "rao",
        N:int = 8,
        expansion_radius_ratio = 0.382,
        clearance:float = 0.1*10**-3,
        **kwargs
    ):
        self.fuel = fuel
        self.oxidizer = oxidizer
        self.chamber_pressure = Pc
        self.exit_pressure_1D = Pe
        self.thrust = thrust
        self.ambient_pressure = Pa
        self.contour = contour.casefold()
        self.xpoints = np.ndarray([])
        self.ypoints = np.ndarray([])
        self.divisions = N
        self.nozzle_OD = nozzle_OD
        self.nozzle_OR = nozzle_OD/2
        self.postcomb_LD = postcomb_LD
        self.angle_converging = angle_converging
        self.expansion_radius_ratio = expansion_radius_ratio
        self.clearance = clearance
        
        self.calculate()
        # Checking for diverging contour type
        if self.contour == "rao":
            self.rao()
        else:
            print("Invalid Contour")
        
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
        
        if dimensions is True:
            self.generate_converging()
            self.model()
            
    
    # # Determine area ratios
    # def area_ratio(self, mach):
    #     k = self.ratio_of_heats
    #     return (
    #         1 / mach * ((2 + (k - 1) * mach**2) / (k + 1)) ** ((k + 1) / (2 * k - 2))
    #     )
        
    # def area_ratio_derivative(self, mach):
    #     k = self.ratio_of_heats
    #     return (
    #     (2 * mach**2 - 2)
    #     * (((k - 1) * mach**2 + 2) / (k + 1)) ** ((k + 1) / (2 * k - 2))
    #     / ((k - 1) * mach**4 + 2 * mach**2)
    #     )
        
    # # Determine mach from area ratios
    # def mach_from_area(self, area, is_subsonic: bool = False):
    #     # General calculation variables
    #     max_iterations = 32
    #     threshold = 10**-16
    #     area_ratio = area / self.area_throat
    #     if is_subsonic:
    #         mach_guess = 1 / area_ratio  # This just seems to work
    #         for i in range(max_iterations):
    #             error = self.area_ratio(mach_guess) - area_ratio
    #             if all(abs(error)) < threshold:
    #                 break
    #             else:
    #                 slope = self.area_ratio_derivative(mach_guess)
    #                 mach_guess -= error / slope
    #     else:
    #         # Uses a 2nd order Taylor series at M=1 to calculate initial guesses for supersonic mach at area ratio
    #         a = 2 / (self.ratio_specific_heat + 1)
    #         mach_guess = 1 + sqrt(abs(area_ratio - 1)) / sqrt(a)
    #         for i in range(max_iterations):
    #             error = self.area_ratio(mach_guess) - area_ratio
    #             if all(abs(error)) < threshold:
    #                 break
    #             else:
    #                 slope = self.area_ratio_derivative(mach_guess)
    #                 mach_guess -= error / slope

    #     return mach_guess
    
    
    def calculate(self):
        # 1D denotes values calculated using the 1D equations, which will be different to those using higher dimensional approahes such as method of
        # characteristics.
        
        # Adding new fuels
        for key in cea_cards.keys():
            # print(key)
            # print(cea_cards[key])
            add_new_fuel(key, cea_cards[key])
        
        # Creation of CEA object
        self.prop = CEA_Obj(oxName=self.oxidizer, fuelName=self.fuel)
        
        # Unit conversion
        self.Pc_psi = self.chamber_pressure*0.0001450377
        self.Pa_psi = self.ambient_pressure*0.0001450377
        
        # Stoichiometri mixture ratio
        self.mixture_ratio = self.prop.getMRforER(ERphi=True)
        
        # self.throat_pressure = self.chamber_pressure*((self.ratio_of_heats + 1)/2)**(self.ratio_of_heats/(1 - self.ratio_of_heats))
        
        self.area_ratio_1D = self.prop.get_eps_at_PcOvPe(Pc=self.Pc_psi, MR=self.mixture_ratio, PcOvPe=self.chamber_pressure/self.exit_pressure_1D, frozen=1)
        
        # self.exit_vel_1D = np.sqrt(
        #     2*self.ratio_of_heats/
        #     (self.mol_mass*(self.ratio_of_heats - 1))*R*self.chamber_temp*(1 - (self.exit_pressure_1D/self.chamber_pressure)**((self.ratio_of_heats - 1)/self.ratio_of_heats))
        #     )
        
        # self.throat_mass_flux = (self.chamber_pressure*self.ratio_of_heats*np.sqrt(
        #     (2/(self.ratio_of_heats + 1))**((self.ratio_of_heats + 1)/(self.ratio_of_heats - 1))
        # )/np.sqrt(self.ratio_of_heats*R*self.chamber_temp/self.mol_mass))
        
        # self.isp_m_s = (self.exit_vel_1D*self.throat_mass_flux + self.area_ratio_1D*(self.exit_pressure_1D - self.ambient_pressure))/self.throat_mass_flux
        
        # self.isp_s = self.isp_m_s/g
        
        self.performance = self.prop.estimate_Ambient_Isp(self.Pc_psi, self.mixture_ratio, self.area_ratio_1D, self.Pa_psi, frozen=1)
        self.isp_s = self.performance[0]
        self.isp_m_s = self.isp_s*g
        self.mass_flow = self.thrust/self.isp_m_s
        self.cstar = self.prop.get_Cstar(self.Pc_psi, self.mixture_ratio)*0.3048
        self.throat_mass_flux = self.cstar*self.prop.get_Densities(self.Pc_psi, self.mixture_ratio, self.area_ratio_1D, frozen=1)[1]*16.018463
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
        print(f"Exit Condition (Psia): {self.performance[1]}")
        print("")
        
    def model(self, export:bool=False, resolution=128):
        # Generating points for CAD
        model_xpoints = [(self.xpoints[0] - self.sheath_length)*1000] + [(self.xpoints[0] - self.sheath_length)*1000] + [self.xpoints[0]*1000] + (self.xpoints*1000).tolist() + [self.xpoints[-1]*1000] + [(self.xpoints[-1] - self.plate_thickness)*1000] + [(self.xpoints[-1] - self.plate_thickness)*1000] + [self.xpoints[0]*1000]
        model_ypoints = [self.nozzle_OR*1000] + [(self.grain.outer_radius + self.clearance)*1000] + [(self.grain.outer_radius + self.clearance)*1000] + (self.ypoints*1000).tolist() + [(self.nozzle_OR - self.lip_thickness)*1000] + [(self.nozzle_OR - self.lip_thickness)*1000] + [self.nozzle_OR*1000] + [self.nozzle_OR*1000]
        model_pts = []
        
        # Generating points for sim
        self.x_values = np.linspace(0, self.xpoints[-1], resolution) # X values of grid points
        self.y_values = np.linspace(0, self.nozzle_OD, resolution) # Y values of grid points
        self.gas_side_y_values = np.interp(self.x_values, self.xpoints, self.ypoints) # Y values of bottom edge
        self.grid_x, self.grid_y = np.meshgrid(self.x_values, self.y_values) # Generates 2d arrays of x values and y values of the grid points
        self.grid_wall_rows = np.zeros(resolution, dtype=int)
        self.grid_wall_columns = np.arange(resolution, dtype=int)
        
        # Clears points that lie outside nozzle, get indicies of wall
        for column in range(resolution):
            external_indices = np.nonzero(self.grid_y[:, column] < self.gas_side_y_values[column])
            self.grid_y[external_indices, column] = NaN
            self.grid_x[external_indices, column] = NaN
            
            # get indicies of points ON GRID that lie on the wall
            self.grid_wall_rows[column] = np.min(np.argwhere(self.grid_y[:, column] > self.gas_side_y_values[column]))
        # Snapping y values to contour
        self.grid_y[self.grid_wall_rows, self.grid_wall_columns] = self.gas_side_y_values
        # get indicies of CELLS that lie on the wall
        self.cell_wall_rows = self.grid_wall_rows[:-1]
        self.cell_wall_columns = self.grid_wall_columns[:-1]
        
        # Generating model
        for i in range(len(model_xpoints)):
            model_pts.append((model_xpoints[i], 0, model_ypoints[i]))
        self.geometry = cq.Workplane("front").polyline(model_pts).close().revolve(axisStart=(0,0,0), axisEnd=(1,0,0), angleDegrees=360)
        # Exporting to step
        if export:
            cq.exporters.export(self.geometry, "nozzle.step")
        
    def rao(self):
        # temp rao nozzle code
        
        # define some required variables
        length_fraction = 1
        parabola_angle_initial = 18
        parabola_angle_final = 9
        angle_div = 15
        
        # determine length of nozzle
        self.length = length_fraction * (
            ((sqrt(self.area_ratio_1D) - 1) * self.throat_radius)
            / (tan(radians(angle_div)))
        )
        
        # converging throat section
        # t1  = linspace(radians(-135), radians(-90), self.divisions, False)
        
        # x_conv = 1.5 * self.throat_radius * cos(t1)
        # y_conv = (1.5 * self.throat_radius * sin(t1)) + (2.5 * self.throat_radius)
        
        # diverging throat section
        t2 = linspace(radians(-90), radians(parabola_angle_initial - 90), self.divisions, False)
        
        x_div = self.expansion_radius_ratio * self.throat_radius * cos(t2)
        y_div = (self.expansion_radius_ratio * self.throat_radius * sin(t2)) + ((1 + self.expansion_radius_ratio) * self.throat_radius)
        
        # start of parabola
        Nx = (
            self.expansion_radius_ratio
            * self.throat_radius
            * cos(radians(parabola_angle_initial) - radians(90))
        )
        
        Ny = (
            self.expansion_radius_ratio
            * self.throat_radius
            * sin(radians(parabola_angle_initial) - radians(90))
        ) + ((1 + self.expansion_radius_ratio) * self.throat_radius)   
        
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
            + ((t3**2) * self.exit_radius_1D)
        )
          
        # add contour to x/y points
        self.xpoints = concatenate((x_div, x_para))
        self.ypoints = concatenate((y_div, y_para))
        
    def generate_converging(self):
        e = (self.grain.outer_radius - self.lip_thickness - self.throat_radius)/np.tan(self.angle_converging)
        
        d = 3/2*self.expansion_radius_ratio*self.throat_radius*np.tan(self.angle_converging)
        b = d
        c = e - b/2 - 5/8*d
        self.radius_postcomb = self.grain.outer_radius - self.lip_thickness
        self.length_postcomb = self.postcomb_LD*(self.grain.outer_diameter - 2*self.lip_thickness)
        L = self.length_postcomb + b + c + d
        x_postccomb = np.linspace(0, self.length_postcomb, self.divisions, endpoint=False)
        x_upstream_quartic = np.linspace(self.length_postcomb, self.length_postcomb + b, self.divisions, endpoint=False)
        x_conical = np.linspace(self.length_postcomb + b, self.length_postcomb + b + c, 1, endpoint=False)
        x_downstream_quartic = np.linspace(self.length_postcomb + b + c, self.length_postcomb + b + c + d, self.divisions, endpoint=False)
        y_chamber = np.full(self.divisions, self.radius_postcomb)
        y_upstream_quartic = self.radius_postcomb - b * tan(
            self.angle_converging
        ) / 2 * ((x_upstream_quartic - self.length_postcomb) / b) ** 3 * (
            2 - (x_upstream_quartic - self.length_postcomb) / b
        )
        y_conical = (
            self.radius_postcomb
            + (self.length_postcomb + b / 2) * tan(self.angle_converging)
            - x_conical * tan(self.angle_converging)
        )
        y_downstream_quartic = (L - x_downstream_quartic) ** 2 / (
            12 * self.expansion_radius_ratio * self.throat_radius
        ) * (6 - ((L - x_downstream_quartic) / d) ** 2) + self.throat_radius
        x_converging = np.concatenate(
            (x_postccomb, x_upstream_quartic, x_conical, x_downstream_quartic)
        )
        y_converging = np.concatenate(
            (y_chamber, y_upstream_quartic, y_conical, y_downstream_quartic)
        )
        self.xpoints = np.concatenate((x_converging, self.xpoints + L))
        self.ypoints = np.concatenate((y_converging, self.ypoints))
        old_xpoints = self.xpoints
        self.length = self.xpoints[-1]
        
    def plot(self):
        plt.plot(self.xpoints, self.ypoints)
        plt.axis('equal')
        plt.show()
        
        
#nozzle = Nozzle(Pc=2*10**6, Tc=3000, thrust=300, M=0.026, mix_ratio=8, y=1.2, ox_den=700)

if __name__ == "__main__":
    nozz = Nozzle(Pc=3*10**6, Tc=3391.91, thrust=300, M=0.026041, mix_ratio=5.3, y=1.2593)
    nozz.plot()