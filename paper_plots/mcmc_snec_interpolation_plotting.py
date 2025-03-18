# import emcee
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import paper_plots.snec_model_interpolator_plotting_sondos as interp
import pandas as pd
import os

# multiplicative factor for SNEC's predicted photospheric velocities to Fe II velocities,
# found to be equal about 1.4, by fitting to all the SNe in our sample
T_thresh = 10 ** 3.75
models = {}


def initialize_empty_models_dict(ranges_dict):
    for data_type in ['lum', 'veloc', 'mag']:
        models[data_type] = {}
        for Mzams in ranges_dict['Mzams']:
            models[data_type][Mzams] = {}
            for Ni in ranges_dict['Ni']:
                models[data_type][Mzams][Ni] = {}
                for E in ranges_dict['E']:
                    models[data_type][Mzams][Ni][E] = {}
                    for R in ranges_dict['R']:
                        models[data_type][Mzams][Ni][E][R] = {}
                        for K in ranges_dict['K']:
                            models[data_type][Mzams][Ni][E][R][K] = {}
                            for Mix in ranges_dict['Mix']:
                                models[data_type][Mzams][Ni][E][R][K][Mix] = None


def get_surrouding_values(requested, ranges_dict):
    params = list(ranges_dict.keys())
    surrouding_values = {param: [] for param in params}
    for i in range(len(requested)):
        param_range = np.array(ranges_dict[params[i]])
        below = np.max(param_range[param_range <= requested[i]])
        above = np.min(param_range[param_range >= requested[i]])
        surrouding_values[params[i]] = [below, above]
    return surrouding_values


def load_model(Mzams, Ni, E, R, K, Mix, data_type, extend_tail=False):
    if R == 0 or K == 0:
        R = 0
        K = 0
    name = 'M' + str(Mzams) + \
           '_Ni' + str(Ni) + \
           '_E' + str(E) + \
           '_Mix' + str(Mix) + \
           '_R' + str(R) + \
           '_K' + str(K)
    if data_type == 'lum':
        modelpath = os.path.join('..', 'SNEC_models', name, 'lum_observed.dat')
        if os.stat(modelpath).st_size < 10 ** 5:
            return 'failed SN'
        else:
            snec_model = pd.read_csv(modelpath)
            time_col = snec_model['t_from_discovery'] / 86400  # sec to days
            interp_days = np.linspace(0, 200, 2001)
            snec_model = np.interp(interp_days, time_col, snec_model['Lum'])
            if extend_tail is not False:
                last_30d_x = interp_days[-100:]
                last_30d_y = snec_model[-100:]
                last_30d_ylog = np.log(last_30d_y)
                tail_poly1d = np.poly1d(np.polyfit(last_30d_x, last_30d_ylog, deg=1))
                extension_days = np.linspace(200.1, 200+extend_tail, int(10*extend_tail))
                extension_lumlog = np.array([tail_poly1d(extension_days[i]) for i in range(len(extension_days))])
                extension_lum = np.exp(extension_lumlog)
                snec_model = np.concatenate((snec_model, extension_lum))
            return snec_model
    elif data_type == 'veloc':
        modelpath = os.path.join('../../../..', 'SNEC_models', name, 'vel_Fe.dat')
        if os.stat(modelpath).st_size < 10 ** 4:
            return 'failed SN'
        else:
            snec_model = pd.read_csv(modelpath)
            interp_days = np.linspace(0, 200, 2001)
            snec_model = np.interp(interp_days, snec_model['t_from_discovery'],
                                   snec_model['vel'])
            return snec_model
    elif data_type == 'mag':
        modelpath = os.path.join('../../../..', 'SNEC_models', name, 'magnitudes_pys.dat')
        lumpath = os.path.join('../../../..', 'SNEC_models', name, 'lum_observed.dat')
        if os.stat(lumpath).st_size < 10 ** 5:
            print('failed SN')
            return 'failed SN'
        else:
            mag_file = pd.read_csv(modelpath, names = ['time', 'Teff', 'PTF_R_AB', 'u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I'])
            mag_file = mag_file.abs()
            time_col = mag_file['time'] / 86400  # sec to days
            interp_days = np.linspace(0, 200, 2001)
            snec_model_dict = {}
            for filter in ['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']:
                snec_model_dict[filter] = np.interp(interp_days, time_col, mag_file[filter])
            snec_model_dict['time'] = interp_days
            snec_model = pd.DataFrame(snec_model_dict)
            snec_model = snec_model.sort_values('time')
            return snec_model

def load_surrounding_models(requested, ranges_dict, fitting_type, extend_tail=False):
    surrouding_values = get_surrouding_values(requested, ranges_dict)
    for Mzams in surrouding_values['Mzams']:
        for Ni in surrouding_values['Ni']:
            for E in surrouding_values['E']:
                for R in surrouding_values['R']:
                    for K in surrouding_values['K']:
                        for Mix in surrouding_values['Mix']:
                            if 'lum' in fitting_type:
                                if models['lum'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['lum'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'lum', extend_tail)
                            if 'veloc' in fitting_type:
                                if models['veloc'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['veloc'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'veloc')
                            if 'mag' in fitting_type:
                                if models['mag'][Mzams][Ni][E][R][K][Mix] is None:
                                    models['mag'][Mzams][Ni][E][R][K][Mix] = load_model(Mzams, Ni, E, R, K, Mix, 'mag')


def theta_in_range(theta, ranges_dict):
    ranges_list = dict_to_list(ranges_dict)
    truth_value = True
    for i in range(len(ranges_list)):
        truth_value = truth_value and\
                      (ranges_list[i][0] <= theta[i] <= ranges_list[i][-1])
    return truth_value


def dict_to_list(dict):
    l = []
    for key in dict.keys():
        l.append(dict[key])
    return l


def rounded_str(x):
    if not np.isinf(x) and np.abs(x) > 0.0000001:
        rounded = round(x, 2-int(np.floor(np.log10(abs(x)))))
        if rounded > 100:
            rounded = int(rounded)
    else:
        rounded = 0
    return str(rounded)


def plot_lum_interpolator_models(requested, ranges_dict, extend_tail=False):
    initialize_empty_models_dict(ranges_dict)
    f_fit, ax = plt.subplots(figsize=(10, 8))
    ax.axvspan(-2, 30, alpha=0.1, color='grey')
    log_likeli = []
    if extend_tail is not False:
        x_plotting = np.linspace(0, 200+extend_tail, int(1+10*200+extend_tail))
    else:
        x_plotting = np.linspace(0, 200, 2001)
    [Mzams, Ni, E, R, K, Mix, S, T] = requested
    if theta_in_range(requested, ranges_dict) or R == 0:
        x_moved = x_plotting - T
        surrounding_values = get_surrouding_values(requested[0:6], ranges_dict)
        load_surrounding_models(requested[0:6], ranges_dict, 'lum', extend_tail)
        y_fit = interp.snec_interpolator(requested[0:6], surrounding_values, models['lum'], x_moved, extend_tail)

Mzams_range = [9.0, 10.0, 11.0, 13.0, 15.0, 17.0]  # Progenitor ZAMS mass (solar mass)
Ni_range = [0.001, 0.02, 0.07, 0.12]  # Ni mass (solar mass)
E_final_range = [0.1, 0.3, 0.5, 0.9, 1.3, 1.7]  # Explosion energy (foe=10^51 erg)
Mix_range = [2.0, 8.0]  # Ni mixing (solar mass)
R_range = [0, 500, 1000]
K_range = [0, 10, 30, 60]
T_range = [-5, 2]  # time shift from assumed day of explosion (days)
parameter_ranges = {'Mzams': Mzams_range, 'Ni': Ni_range, 'E': E_final_range,
                    'R': R_range, 'K': K_range, 'Mix': Mix_range, 'S': [0.1, 2], 'T': T_range}

# change these values to determined the values of the requested the interpolated model example
Mzams = 14.1       #12.2
Ni =  0.05          #0.01
E =  1.6          #1.11
R =  828.2          #650
K =  31.3          #0.01
Mix =  3        #4.2
S = 1           #1
T = 0
requested_theta = [Mzams, Ni, E, R, K, Mix, S, T]

plot_lum_interpolator_models(requested_theta, parameter_ranges)




'''

def plot_mag_interpolator_models(requested, ranges_dict):
    initialize_empty_models_dict(ranges_dict)
    f_fit, ax = plt.subplots(figsize=(10, 8))
    ax.axvspan(-2, 30, alpha=0.1, color='grey')
    log_likeli = []
    x_plotting = np.linspace(0, 200, 2001)
    [Mzams, Ni, E, R, K, Mix, S, T] = requested
    if theta_in_range(requested, ranges_dict) or R == 0:
        load_surrounding_models(requested[0:6], ranges_dict, 'mag')
        x_moved_all = x_plotting - T
        surrounding_values = get_surrouding_values(requested[0:6], ranges_dict)
        y_fit = interp.snec_interpolator(requested[0:6], surrounding_values, models['mag'], x_moved_all)



def plot_veloc_interpolator_models(requested, ranges_dict):
    initialize_empty_models_dict(ranges_dict)
    f_fit, ax = plt.subplots(figsize=(10, 8))
    log_likeli = []
    [Mzams, Ni, E, R, K, Mix, S, T] = requested
    if theta_in_range(requested, ranges_dict) or R == 0:
        load_surrounding_models(requested[0:6], ranges_dict, 'veloc')
        surrounding_values = get_surrouding_values(requested[0:6], ranges_dict)
        x_plotting = np.linspace(0, 200, 2001)
        x_moved = x_plotting - T
        y_fit = interp.snec_interpolator(requested[0:6], surrounding_values, models['veloc'], x_moved)


'''