from constants import *
import numpy as np
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
        cd:float = 1,
        Pe:float = P_sl
    ):
        self.chamber_pressure = Pc
        self.exit_pressure_1D = Pe
        self.chamber_temp = Tc
        self.mixture_ratio = mix_ratio
        self.ratio_of_heats = y
        self.mol_mass = M
        self.thrust = thrust
        self.ox_density = ox_den
        self.discharge_coeff = cd
        
        # calcs
        self.pressure_ratio = self.chamber_pressure/self.exit_pressure_1D
        
        self.area_ratio_1D = 1/(((self.ratio_of_heats + 1)/2)**(1/(self.ratio_of_heats - 1))*(self.exit_pressure_1D/self.chamber_pressure)**(1/self.ratio_of_heats)
                              *np.sqrt((self.ratio_of_heats + 1)/(self.ratio_of_heats - 1)*(1 - (self.exit_pressure_1D/self.chamber_pressure)**((self.ratio_of_heats - 1)/self.ratio_of_heats)))
                              )
        
        self.exit_vel_1D = np.sqrt(
            2*self.ratio_of_heats/
            (self.mol_mass*(self.ratio_of_heats - 1))*R*self.chamber_temp*(1 - (self.exit_pressure_1D/self.chamber_pressure)**((self.ratio_of_heats - 1)/self.ratio_of_heats))
            )
        #self.isp_m_s = (self.exit_vel_1D*self.mass_flux())
        
        self.ox_flow = self.thrust/self.exit_vel_1D
        
    # funcs
    def mass_flux(self, pressure):
        return self.chamber_pressure*self.ratio_of_heats*np.sqrt(
            (2/(self.ratio_of_heats + 1))**((self.ratio_of_heats + 1)/(self.ratio_of_heats - 1))
        )/np.sqrt(self.ratio_of_heats*R*self.chamber_temp)

engine = Engine(Pc=2*10**6, Tc=3000, thrust=300, M=0.026, mix_ratio=8, y=1.2)
print(engine.ox_flow)