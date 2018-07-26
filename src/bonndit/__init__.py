# -*- coding: utf-8 -*-

"""Top-level package for bonndit."""

__author__ = """Olivier Morelle"""
__email__ = 'morelle@uni-bonn.de'
__version__ = '0.1.0'

from dipy.reconst.dti import TensorModel

from .cylkurtosis import CylKurtosisModel
from .io import load
from .shoremt import ShoreFitMt, ShoreModelMt
from .shorest import ShoreModel

# In future we want to add all models from dipy and bonndit
dwi_models = {  # "mtshore": mtShoreModel,
    "cylkurtosis": CylKurtosisModel,
    "tensor": TensorModel}

conv_frameworks = {"mtshore": ShoreFitMt, }
