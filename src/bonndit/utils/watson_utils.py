import numpy as np
import os

from bonndit.directions.fodfapprox import approx_all_spherical
from bonndit.utils.esh import sym_to_esh, esh_to_sym
from bonndit.utils.storage import nrrd_to_meta3d
from bonndit.utils.fields import save_tensor

EPS_L = 0.01
CONF_90 = lambda x : 9.22022495e+01 + x*(-1.19665463e+01 + x*(8.57242880e-01 + x*(-3.50226232e-02 + x*(8.84091225e-04 + x*(-1.44771428e-05 + x*(1.57390844e-07 + x*(-1.13575538e-09 + x*(5.27740057e-12 + x*(-1.44255882e-14 + x*(1.68532222e-17 + x*(1.13446824e-20 + x*(-3.69532677e-23))))))))))))
O_LIMIT = np.array([0,0,0,0,0.02,0,0.04,0,0.08])

def fodf_to_watson_sh(fodf_data, shorder):
    fodf_track_data = np.moveaxis(fodf_data, 0, -1)[:,:,:,1:]

    # convert to sh
    shm_coeff = np.zeros_like(fodf_track_data)
    for i in range(fodf_track_data.shape[0]):
        for j in range(fodf_track_data.shape[1]):
            for k in range(fodf_track_data.shape[2]):
                shm_coeff[i,j,k] = sym_to_esh(fodf_track_data[i,j,k])

    # reorder sh coeffs
    flipped_shm_coeff = v_flip_sh_signal(shm_coeff.reshape(-1, shm_coeff.shape[3]), shorder)
    shm_coeff = flipped_shm_coeff.reshape(shm_coeff.shape)
    return shm_coeff

def watson_sh_to_fodf(sh_data, shorder):
    # reorder sh coeffs
    flipped_sh_data = v_flip_sh_signal(sh_data.reshape(-1, sh_data.shape[3]), shorder)
    sh_data = flipped_sh_data.reshape(sh_data.shape)

    # convert to fodf
    fodf_coeff = np.zeros_like(sh_data)
    for i in range(sh_data.shape[0]):
        for j in range(sh_data.shape[1]):
            for k in range(sh_data.shape[2]):
                fodf_coeff[i,j,k] = esh_to_sym(sh_data[i,j,k])

    return fodf_coeff

def lowrankapprox_from_fodf(fodf_data, rank, verbose=False):
    data = fodf_data.reshape((fodf_data.shape[0], -1))
    print(data.shape[1])
    output = np.zeros((4, rank, data.shape[1]))
    approx_all_spherical(output, np.float64(fodf_data), int(0), np.float32(0), rank, verbose)
    output = output.reshape((4, rank) + fodf_data.shape[1:])
    return nrrd_peaks_split(output)
    
def nrrd_peaks_split(nrrd_peaks):
    nrrd_peaks = np.moveaxis(nrrd_peaks, [0, 1], [-1, -2])
    peak_dirs = nrrd_peaks[:,:,:,:,1:]
    peak_values = nrrd_peaks[:,:,:,:,0]
    return peak_values, peak_dirs

def create_missing_folders(path_with_file):
    folders = os.path.split(path_with_file)[0]
    if folders != '' and not os.path.exists(folders):
        os.makedirs(folders)

def save_as_fodf_nrrd(filename, fodf, header):
    meta = nrrd_to_meta3d(header)
    save_tensor(filename, fodf, meta=meta)


def v_flip_sh_signal(sh, shorder):
    if shorder == 4:
        return v_flip_sh_signal_o4(sh)
    elif shorder == 6:
        return v_flip_sh_signal_o6(sh)
    elif shorder == 8:
        return v_flip_sh_signal_o8(sh)

def v_flip_sh_signal_o4(sh):
    mapped_sh = sh.copy()
    mapped_sh[:,5] = sh[:,1]
    mapped_sh[:,4] = sh[:,2]
    mapped_sh[:,2] = sh[:,4]
    mapped_sh[:,1] = sh[:,5]
    mapped_sh[:,14] = sh[:,6]
    mapped_sh[:,13] = sh[:,7]
    mapped_sh[:,12] = sh[:,8]
    mapped_sh[:,11] = sh[:,9]
    mapped_sh[:,9] = sh[:,11]
    mapped_sh[:,8] = sh[:,12]
    mapped_sh[:,7] = sh[:,13]
    mapped_sh[:,6] = sh[:,14]

    return mapped_sh

def flip_sh_signal_o4(sh):
    mapped_sh = sh.copy()
    mapped_sh[5] = sh[1]
    mapped_sh[4] = sh[2]
    mapped_sh[2] = sh[4]
    mapped_sh[1] = sh[5]
    mapped_sh[14] = sh[6]
    mapped_sh[13] = sh[7]
    mapped_sh[12] = sh[8]
    mapped_sh[11] = sh[9]
    mapped_sh[9] = sh[11]
    mapped_sh[8] = sh[12]
    mapped_sh[7] = sh[13]
    mapped_sh[6] = sh[14]

    return mapped_sh

def v_flip_sh_signal_o6(sh):
    mapped_sh = sh.copy()
    mapped_sh[:,5]  = sh[:,1]
    mapped_sh[:,4]  = sh[:,2]
    mapped_sh[:,2]  = sh[:,4]
    mapped_sh[:,1]  = sh[:,5]
    mapped_sh[:,14] = sh[:,6]
    mapped_sh[:,13] = sh[:,7]
    mapped_sh[:,12] = sh[:,8]
    mapped_sh[:,11] = sh[:,9]
    mapped_sh[:,9]  = sh[:,11]
    mapped_sh[:,8]  = sh[:,12]
    mapped_sh[:,7]  = sh[:,13]
    mapped_sh[:,6]  = sh[:,14]
    mapped_sh[:,20] = sh[:,15]
    mapped_sh[:,19] = sh[:,16]
    mapped_sh[:,18] = sh[:,17]
    mapped_sh[:,17] = sh[:,18]
    mapped_sh[:,16] = sh[:,19]
    mapped_sh[:,15] = sh[:,20]
    mapped_sh[:,27] = sh[:,22]
    mapped_sh[:,26] = sh[:,23]
    mapped_sh[:,25] = sh[:,24]
    mapped_sh[:,24] = sh[:,25]
    mapped_sh[:,23] = sh[:,26]
    mapped_sh[:,22] = sh[:,27]

    return mapped_sh

def flip_sh_signal_o6(sh):
    mapped_sh = sh.copy()
    mapped_sh[5]  = sh[1]
    mapped_sh[4]  = sh[2]
    mapped_sh[2]  = sh[4]
    mapped_sh[1]  = sh[5]
    mapped_sh[14] = sh[6]
    mapped_sh[13] = sh[7]
    mapped_sh[12] = sh[8]
    mapped_sh[11] = sh[9]
    mapped_sh[9]  = sh[11]
    mapped_sh[8]  = sh[12]
    mapped_sh[7]  = sh[13]
    mapped_sh[6]  = sh[14]
    mapped_sh[20] = sh[15]
    mapped_sh[19] = sh[16]
    mapped_sh[18] = sh[17]
    mapped_sh[17] = sh[18]
    mapped_sh[16] = sh[19]
    mapped_sh[15] = sh[20]
    mapped_sh[27] = sh[22]
    mapped_sh[26] = sh[23]
    mapped_sh[25] = sh[24]
    mapped_sh[24] = sh[25]
    mapped_sh[23] = sh[26]
    mapped_sh[22] = sh[27]

    return mapped_sh

def v_flip_sh_signal_o8(sh):
    mapped_sh = sh.copy()
    mapped_sh[:,5]  = sh[:,1]
    mapped_sh[:,4]  = sh[:,2]
    mapped_sh[:,2]  = sh[:,4]
    mapped_sh[:,1]  = sh[:,5]
    mapped_sh[:,14] = sh[:,6]
    mapped_sh[:,13] = sh[:,7]
    mapped_sh[:,12] = sh[:,8]
    mapped_sh[:,11] = sh[:,9]
    mapped_sh[:,9]  = sh[:,11]
    mapped_sh[:,8]  = sh[:,12]
    mapped_sh[:,7]  = sh[:,13]
    mapped_sh[:,6]  = sh[:,14]
    mapped_sh[:,20] = sh[:,15]
    mapped_sh[:,19] = sh[:,16]
    mapped_sh[:,18] = sh[:,17]
    mapped_sh[:,17] = sh[:,18]
    mapped_sh[:,16] = sh[:,19]
    mapped_sh[:,15] = sh[:,20]
    mapped_sh[:,27] = sh[:,22]
    mapped_sh[:,26] = sh[:,23]
    mapped_sh[:,25] = sh[:,24]
    mapped_sh[:,24] = sh[:,25]
    mapped_sh[:,23] = sh[:,26]
    mapped_sh[:,22] = sh[:,27]
    mapped_sh[:,35] = sh[:,28]
    mapped_sh[:,34] = sh[:,29]
    mapped_sh[:,33] = sh[:,30]
    mapped_sh[:,32] = sh[:,31]
    mapped_sh[:,31] = sh[:,32]
    mapped_sh[:,30] = sh[:,33]
    mapped_sh[:,29] = sh[:,34]
    mapped_sh[:,28] = sh[:,35]
    mapped_sh[:,44] = sh[:,37]
    mapped_sh[:,43] = sh[:,38]
    mapped_sh[:,42] = sh[:,39]
    mapped_sh[:,41] = sh[:,40]
    mapped_sh[:,40] = sh[:,41]
    mapped_sh[:,39] = sh[:,42]
    mapped_sh[:,38] = sh[:,43]
    mapped_sh[:,37] = sh[:,44]

    return mapped_sh