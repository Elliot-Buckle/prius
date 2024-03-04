from engine import Engine
class Nozzle:
    class Contours:
        def __init__(
            self,
            engine:Engine,
            contour:str = 'moc'
        ):
            self.engine = engine
            self.contour = contour
            
            #Calculations