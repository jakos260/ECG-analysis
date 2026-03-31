from carto3reader.CartoMapDataPaths import CartoMapDataPaths
from carto3reader.CartoMapModel import CartoMapModel

class CartoMap(CartoMapDataPaths, CartoMapModel):
    def __init__(self, data_path, measurement_name):
        CartoMapDataPaths.__init__(self, data_path, measurement_name)
        CartoMapModel.__init__(self)
        self._load_master()


class CartoMaps(CartoMap):
    pass