import os
from os.path import join


class Paths:
    def __init__(self):
        self.sw_install = '/data/Evo/Research/ZoomQA'
        self.model_path = join(self.sw_install, "model/")
        model_name = os.listdir(self.model_path)[0]
        self.model_path = join(self.model_path, model_name)


PATHS = Paths()
