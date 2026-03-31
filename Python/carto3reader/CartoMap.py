from carto3reader.CartoMapDataPaths import CartoMapDataPaths
from carto3reader.CartoMapModel import CartoMapModel
from carto3reader.EcgReader import EcgReader

class CartoMap(CartoMapDataPaths, CartoMapModel):
    def __init__(self, data_path, measurement_name):
        CartoMapDataPaths.__init__(self, data_path, measurement_name)
        CartoMapModel.__init__(self)
        self.ecg_reader = EcgReader(self)
        self._load_master()


class CartoMaps(CartoMap):
    pass