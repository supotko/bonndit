import os

dir_path = os.path.dirname(os.path.realpath(__file__))


DECONVOLUTION_DIR = os.path.join(dir_path, 'data/deconvolution')
DECONVOLUTION_RESULTS_DIR = os.path.join(DECONVOLUTION_DIR, 'results')

DWMRI_DUMMY_DATA = os.path.join(DECONVOLUTION_DIR, 'data.nii.gz')
SHORE_FIT_TEST = os.path.join(DECONVOLUTION_RESULTS_DIR, 'response.npz')
# Response computed with michi-temp shore-response (rotated to worldspace)
SHORE_FIT_PRECOMPUTED = os.path.join(DECONVOLUTION_DIR, 'response_complete.npz')

ODF_RESULT_HPSD = os.path.join(DECONVOLUTION_RESULTS_DIR, 'odfhpsd.nrrd')
ODF_RESULT_NONNEG = os.path.join(DECONVOLUTION_RESULTS_DIR, 'odfnonneg.nrrd')
ODF_RESULT_NO_CONSTRAINT = os.path.join(DECONVOLUTION_RESULTS_DIR, 'odfnone.nrrd')
