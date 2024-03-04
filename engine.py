from constants import *
from injector import Injector
from nozzle import Nozzle
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
        Pe:float = P_sl,
        Pa:float = P_sl
    ):
        self.nozzle = Nozzle(Pc, Tc, thrust, M, mix_ratio, y, ox_den, cd, Pe, Pa)
        self.injector = Injector(self.nozzle)

engine = Engine(Pc=2*10**6, Tc=3000, thrust=300, M=0.026, mix_ratio=8, y=1.2, ox_den=700)
#print(engine.nozzle.isp_s)