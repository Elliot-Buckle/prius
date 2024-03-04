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
        Pe:float = P_sl
    ):
        self.chamber_pressure = Pc
        self.exit_pressure_1D = Pe
        self.chamber_temp = Tc
        self.mixture_ratio = mix_ratio
        self.ratio_of_heats = y
        self.mol_mass = M
        self.thrust = thrust
        
        #calcs
        self.exit_vel_1D = np.sqrt(
            2*self.ratio_of_heats/(self.mol_mass*(self.ratio_of_heats - 1))*R*self.chamber_temp*(1 - (self.exit_pressure_1D/self.chamber_pressure)**((self.ratio_of_heats - 1)/self.ratio_of_heats))
            )
        self.ox_flow = self.thrust/self.exit_vel_1D

engine = Engine(Pc=2*10**6, Tc=3000, thrust=300, M=0.026, mix_ratio=8, y=1.2)
print(engine.ox_flow)