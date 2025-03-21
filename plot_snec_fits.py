import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import os
import mcmc_snec
import re
import copy
import colorful_corner_plot as color_corner

T_thresh = 10 ** 3.5   #(recombination)
extend_tail = False
# filters = ['u', 'g', 'r', 'i', 'z', 'U', 'B', 'V', 'R', 'I']
filters = ['g', 'r', 'i', 'z', 'B', 'V', 'R', 'I']
colors = {'u': 'indigo', 'g': 'green', 'r': 'red', 'i': 'dimgrey', 'z': 'black', 'U': 'purple',
          'B': 'blue', 'V': 'darkviolet', 'R': 'darkorange', 'I': 'cyan'}



def add_likelihood_to_file(model_name, fit_type, likeli, output_dir):
    filepath = os.path.join(output_dir, 'log_likelihoods.txt')
    if not os.path.exists(filepath):
        f = pd.DataFrame({'model_name':[], 'fit_type':[], 'log_likelihood':[]})
        f.to_csv(filepath)
    f = pd.read_csv(filepath, index_col=0)
    f = f.append(pd.DataFrame({'model_name':model_name, 'fit_type':fit_type, 'log_likelihood':likeli},index=[0]),ignore_index=True)
    f.to_csv(filepath)


def import_ranges(list, run_params):
    ranges = []
    for name in list:
        if 'range' in name:
            content = re.sub(r'[\[\]\s]', '', run_params.iloc[0][name]).split(',')
            if 'R_' in name or 'K_' in name:
                content = [int(content[i]) for i in range(len(content))]
            else:
                content = [float(content[i]) for i in range(len(content))]
            ranges.append(content)
    ranges_dict = {'Mzams': ranges[0], 'Ni': ranges[1], 'E': ranges[2],
                   'R': ranges[3], 'K': ranges[4], 'Mix': ranges[5],
                   'S': ranges[6], 'T': ranges[7]}
    return ranges_dict


def get_param_results_dict(sampler_df, ranges_dict):
    params = list(ranges_dict.keys())
    dict = {}
    for param in params:
        avg = np.average(sampler_df[param])
        sigma_lower, sigma_upper = np.percentile(sampler_df[param], [16, 84])
        dict[param] = avg
        dict[param + '_lower'] = avg - sigma_lower
        dict[param + '_upper'] = sigma_upper - avg
    return dict


def rounded_str(x):
    if not np.isinf(x) and np.abs(x) > 0.0000001:
        rounded = round(x, 2-int(np.floor(np.log10(abs(x)))))
        if rounded > 100:
            rounded = int(rounded)
    else:
        rounded = 0
    return str(rounded)


def result_text_from_dict(param_dict, ranges_dict):
    params = list(ranges_dict.keys())
    res_text = ''
    for param in params:
        if (param != 'K') & (param != 'R'):
            res_text += param + ': ' + rounded_str(param_dict[param]) + r'$\pm$ [' +\
                            rounded_str(param_dict[param+'_lower']) + ',' +\
                            rounded_str(param_dict[param+'_upper']) + ']\n'
    if (param_dict['K'] == 0) & (param_dict['R'] == 0):
        res_text += 'no CSM'
    else:
        res_text += 'K' + ': ' + rounded_str(param_dict['K']) + r'$\pm$ [' + \
                    rounded_str(param_dict['K_lower']) + ',' + \
                    rounded_str(param_dict['K_upper']) + ']\n'
        res_text += 'R' + ': ' + rounded_str(param_dict['R']) + r'$\pm$ [' + \
                    rounded_str(param_dict['R_lower']) + ',' + \
                    rounded_str(param_dict['R_upper']) + ']\n'
    return res_text


def open_reshape_3d_array(output_dir, type, step):
    array_2d = np.genfromtxt(os.path.join(output_dir, type+ '_models_step' + str(step) + '.txt'))
    shape_2d = array_2d.shape
    array_3d = array_2d.reshape((shape_2d[0], 2, shape_2d[1]/2))
    return array_3d


def add_model_martinez(SN_name, ranges_dict, x_data, y_data, dy_data, fig_type, ax):
    martinez_path = os.path.join('SN_data', 'martinez_results_table.csv')
    martinez_df = pd.read_csv(martinez_path)
    results_row = martinez_df.loc[martinez_df['SN'] == SN_name]
    results_row = results_row.iloc[0]
    log_likeli = []
    requested = [results_row[a] for a in ['Mzams','Ni','E','R','K','Mix','S','T']]
    print(requested)
    S = requested[6]
    T = requested[7]
    data_x_moved = x_data - T
    x_plotting = np.linspace(-T, 200 - T, 2001)
    # Restrict x_plotting to the range of x_data for velocity
    #x_plotting_veloc = x_plotting[x_plotting <= x_data.max()] (עד הדאטא)
    x_plotting_veloc = x_plotting[x_plotting <=125]

    y_fit_plotting = mcmc_snec.interp_yfit(requested, ranges_dict, fig_type, x_plotting)
    y_fit_on_data_times = mcmc_snec.interp_yfit(requested, ranges_dict, fig_type, data_x_moved)
    if not isinstance(y_fit_plotting, str):
        if fig_type == 'lum':
            y_fit_plotting = y_fit_plotting * S
            y_fit_on_data_times = y_fit_on_data_times * S
            ax.plot(x_plotting, np.log10(y_fit_plotting), alpha=0.8, color='orange')
        elif fig_type == 'veloc':
            #ax.plot(x_plotting, y_fit_plotting, alpha=0.8, color='orange') #default color='purple'
            # Only plot the orange curve up to the end of the data
            y_fit_plotting_veloc = mcmc_snec.interp_yfit(requested, ranges_dict, fig_type, x_plotting_veloc)
            ax.plot(x_plotting_veloc, y_fit_plotting_veloc, alpha=0.8, color='orange')
        else:
            print('fig_type must be lum or veloc')
        log_likeli.append(
            mcmc_snec.calc_likelihood(data_x_moved, y_data, dy_data, y_fit_on_data_times, normalization=False))
    log_likeli = np.mean(log_likeli)
    #param_dict = dict(results_row.iloc[0])     #replaced this with: #sondos (ask noi to make sure)
    param_dict = {column_name: value for column_name, value in results_row.items()}
    results_text = result_text_from_dict(param_dict, ranges_dict)
    handles, labels = plt.gca().get_legend_handles_labels()
    SNemcee_fit_lengend = Line2D([0], [0], label=results_text, color='orange')
    handles.extend([SNemcee_fit_lengend])
    ax.legend(handles=handles, fontsize=12)
    output_dir = os.path.join('figures', 'martinez')
    add_likelihood_to_file(SN_name+'_martinez', 'lum', log_likeli, output_dir)
    return ax



def plot_lum_with_fit(SN_name, data_dict, sampler_df, ranges_dict, n_walkers, ax, normalization, LumThreshold, add_martinez=False):
    data = data_dict['lum']
    data_x = data['t_from_discovery']
    data_y = data['Lum']
    dy0 = data['dLum0']
    dy1 = data['dLum1']
    log_likeli = []
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_df.iloc[i]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if mcmc_snec.theta_in_range(requested, ranges_dict):
            data_x_moved = data_x - T
            data_tempthresh = copy.deepcopy(data)
            if LumThreshold:
                max_x, temp_fit = mcmc_snec.temp_thresh_cutoff(requested[0:6], ranges_dict, data_x_moved)
                data_tempthresh = data_tempthresh.loc[data_x_moved <= max_x]
                data_x_moved = data_tempthresh['t_from_discovery'] - T
                x_plotting = np.linspace(-T, max_x, int(1 + max_x * 10))
            else:
                x_plotting = np.linspace(-T, 200-T, 2001)
            y_fit_plotting = mcmc_snec.interp_yfit(requested, ranges_dict, 'lum', x_plotting)
            y_fit_on_data_times = mcmc_snec.interp_yfit(requested, ranges_dict, 'lum', data_x_moved)
            if not isinstance(y_fit_plotting, str):
                # multiply whole graph by scaling factor
                y_fit_plotting = y_fit_plotting * S
                y_fit_on_data_times = y_fit_on_data_times * S
                ax.plot(x_plotting, np.log10(y_fit_plotting), alpha=0.1, color='purple')
                log_likeli.append(mcmc_snec.calc_likelihood(data_x_moved, data_tempthresh['Lum'], data_tempthresh['dLum0'], y_fit_on_data_times, normalization))
    log_likeli = np.mean(log_likeli)
    if add_martinez:
        ax = add_model_martinez(SN_name, ranges_dict, data_x, data_y, dy0, 'lum', ax)
    data_dy0 = np.log10(data_y + dy0) - np.log10(data_y)
    data_dy1 = np.log10(data_y + dy1) - np.log10(data_y)
    ax.errorbar(data_x, np.log10(data_y), yerr=[data_dy0, data_dy1], marker='o', linestyle='None', color='k')
    param_dict = get_param_results_dict(sampler_df, ranges_dict)
    results_text = result_text_from_dict(param_dict, ranges_dict)
    handles, labels = plt.gca().get_legend_handles_labels()
    SNemcee_fit_lengend = Line2D([0], [0], label=results_text, color='purple')
    handles.extend([SNemcee_fit_lengend])
    ax.legend(handles=handles, fontsize=10)
    ax.set_xlim(-2, 200)   #here
    ax.set_ylim(top=43.5)
    ax.tick_params(axis='both', which='major', labelsize=10)
    return ax, log_likeli


def plot_veloc_with_fit(SN_name, data_dict, sampler_df, ranges_dict, n_walkers, ax, normalization, LumThreshold, add_martinez=False):
    data = data_dict['veloc']
    data_x = data['t_from_discovery']
    data_y = data['veloc']
    data_dy = data['dveloc']
    log_likeli = []
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_df.iloc[i]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if mcmc_snec.theta_in_range(requested, ranges_dict):
            data_x_moved = data_x - T
            data_tempthresh = copy.deepcopy(data)
            # always truncate by temp thresh for veloc
            max_x, temp_fit = mcmc_snec.temp_thresh_cutoff(requested[0:6], ranges_dict, data_x_moved)
            data_tempthresh = data_tempthresh.loc[data_x_moved <= max_x]
            data_x_moved = data_tempthresh['t_from_discovery'] - T
            x_plotting = np.linspace(-T, max_x, int(1 + max_x * 10))
            y_fit_plotting = mcmc_snec.interp_yfit(requested, ranges_dict, 'veloc', x_plotting)
            y_fit_on_data_times = mcmc_snec.interp_yfit(requested, ranges_dict, 'veloc', data_x_moved)
            if not isinstance(y_fit_plotting, str):
                ax.plot(x_plotting, y_fit_plotting, alpha=0.1, color='purple')
                log_likeli.append(mcmc_snec.calc_likelihood(data_x_moved, data_tempthresh['veloc'], data_tempthresh['dveloc'], y_fit_on_data_times, normalization))
    log_likeli = np.mean(log_likeli)
    if add_martinez:
        ax = add_model_martinez(SN_name, ranges_dict, data_x, data_y, data_dy, 'veloc', ax)
    # real observations for the SN
    ax.errorbar(data_x, data_y, yerr=data_dy, marker='o', linestyle='None', color='k')
    ax.set_xlim(-2, 200)
    ax.tick_params(axis='both', which='major', labelsize=10)
    return ax, log_likeli


def plot_mag_with_fit(SN_name, data_dict, sampler_df, ranges_dict, n_walkers, ax, normalization, LumThreshold, add_martinez):
    data = data_dict['mag']
    filters = list(data['filter'].unique())
    data_x = data['t_from_discovery']
    y_fit_plotting = {}
    y_fit_on_data_times = {}
    log_likeli = []

    # Setting a uniform offset for each filter.
    offsets = {filt: 1.0 * i for i, filt in enumerate(filters)}
    for i in range(n_walkers):
        [Mzams, Ni, E, R, K, Mix, S, T] = sampler_df.iloc[i]
        requested = [Mzams, Ni, E, R, K, Mix, S, T]
        if mcmc_snec.theta_in_range(requested, ranges_dict):
            data_x_moved = data_x - T
            data_tempthresh = copy.deepcopy(data)
            # always truncate by temp thresh for mag
            max_x, temp_fit = mcmc_snec.temp_thresh_cutoff(requested[0:6], ranges_dict, data_x_moved)
            data_tempthresh = data_tempthresh.loc[data_x_moved <= max_x]
            x_plotting = np.linspace(-T, max_x, int(1 + max_x * 10))

            for filt in filters:
                offset = offsets[filt]
                data_filt = data_tempthresh.loc[data_tempthresh['filter'] == filt]
                data_x_filt_moved = data_filt['t_from_discovery'] - T
                data_y_filt = data_filt['abs_mag']
                data_dy_filt = data_filt['dmag']
                y_fit_plotting[filt] = mcmc_snec.interp_yfit(requested, ranges_dict, 'mag', x_plotting, filt)
                y_fit_on_data_times[filt] = mcmc_snec.interp_yfit(requested, ranges_dict, 'mag', data_x_filt_moved, filt)
                # multiply whole graph by scaling factor
                if not isinstance(y_fit_plotting[filt], str):
                    y_fit_plotting[filt] = y_fit_plotting[filt] - 2.5*np.log10(S)
                    y_fit_on_data_times[filt] = y_fit_on_data_times[filt] - 2.5*np.log10(S)
                    ax.plot(x_plotting, y_fit_plotting[filt] + offset, color=colors[filt], alpha=0.1)
                    log_likeli.append(mcmc_snec.calc_likelihood(
                        data_x_filt_moved, data_y_filt, data_dy_filt, y_fit_on_data_times[filt], normalization
                    ))

    log_likeli = np.mean(log_likeli)
    # log_likeli = np.sum(log_likeli)
    for filt in filters:
        offset = offsets[filt]
        data_filt = data.loc[data['filter'] == filt]
        data_x = data_filt['t_from_discovery']
        data_y = data_filt['abs_mag'] - 2.5 * np.log10(S) + offset
        data_dy = data_filt['dmag']

        label = f'{filt} + {offset:.1f}' if offset != 0 else filt
        ax.errorbar(data_x, data_y, yerr=data_dy, marker='o', linestyle='None', label=label, color=colors[filt])

    ax.legend()
    ax.set_xlim(-2, 200)
    ax.invert_yaxis()
    ax.set_ylim(-8, -19)
    ax.tick_params(axis='both', which='major', labelsize=10)

    return ax, log_likeli



def range_bounds(ranges_list):
    tup_list = []
    for i in range(len(ranges_list)):
        tup_list.append((np.min(ranges_list[i]), np.max(ranges_list[i])))
    return tup_list


def overlay_corner_plot(result_paths_list, output_dir, name_list, filename):
    ranges_dicts_list = []
    flat_sampler_list = []
    for result_path in result_paths_list:
        run_params_path = os.path.join(result_path, 'run_parameters.csv')
        run_params = pd.read_csv(run_params_path, index_col=0).T
        ranges_dict = import_ranges(run_params.columns.values, run_params)
        params = ranges_dict.keys()
        n_walkers = int(run_params.iloc[0]['n_walkers'])
        burn_in = int(run_params.iloc[0]['burn_in'])
        flat_sampler_path = os.path.join(result_path, 'flat_sampler.csv')
        flat_sampler_noburnin = pd.read_csv(flat_sampler_path,
                                            names=params,
                                            skiprows=(burn_in - 1) * (n_walkers))
        ranges_dicts_list.append(ranges_dict)
        flat_sampler_list.append(flat_sampler_noburnin)
    labels = list(ranges_dicts_list[0].keys())
    corner_range = [1.] * len(labels)
    color_corner.overlaid_corner(
        flat_sampler_list,
        name_list,
        corner_range,
        labels,
        output_dir, filename)
     #f_corner = corner.corner(sampler_chain_flat, labels=labels, range=corner_range)
     #f_corner.savefig(os.path.join(output_dir, 'corner_plot.png'))


def string_to_bool(mystr):
    if mystr == 'False':
        return False
    else:
        return True


def get_each_walker_result(sampler_chain, ranges_dict, step):
    params = list(ranges_dict.keys())
    dict = {}
    for i in range(len(params)):
        last_results = sampler_chain[:, step:, i]
        avg = np.mean(last_results)
        sigma_lower, sigma_upper = np.percentile(last_results, [16, 84])
        dict[params[i]] = avg
        dict[params[i]+'_lower'] = avg - sigma_lower
        dict[params[i] + '_upper'] = sigma_upper - avg

def chain_plots(result_path, output_dir, first_stage_steps=None):
    # (result_paths_list, output_dir, name_list, filename):
    run_params_path = os.path.join(result_path, 'run_parameters.csv')
    run_params = pd.read_csv(run_params_path, index_col=0).T
    ranges_dict = import_ranges(run_params.columns.values, run_params)
    params = ranges_dict.keys()
    n_walkers = int(run_params.iloc[0]['n_walkers'])
    n_steps = int(run_params.iloc[0]['n_steps'])
    burn_in = int(run_params.iloc[0]['burn_in'])
    flat_sampler_path = os.path.join(result_path, 'flat_sampler.csv')
    sampler_array = np.genfromtxt(flat_sampler_path, delimiter=',')
    sampler_chain = sampler_array.reshape((n_walkers, n_steps, len(params)))
    keys = list(ranges_dict.keys())
    for i in range(len(keys)):
        key = keys[i]
        plt.figure()
        plt.plot(sampler_chain[:, :, i].T, color='k', alpha=0.1)
        plt.xlabel('Step Number')
        plt.ylabel(key)
        plt.axvspan(0, burn_in, alpha=0.1, color='grey')
        if first_stage_steps is not None:
            plt.axvline(x=first_stage_steps, color='black')
            plt.axvspan(first_stage_steps, first_stage_steps+burn_in, alpha=0.1, color='grey')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, key+'.png'))
        plt.savefig(os.path.join(output_dir, key + '.pdf'))


def get_args_from_file(result_path, ax, data_type, add_martinez=False):
    run_params_path = os.path.join(result_path, 'run_parameters.csv')
    run_params = pd.read_csv(run_params_path, index_col=0).T
    ranges_dict = import_ranges(run_params.columns.values, run_params)
    params = ranges_dict.keys()
    SN_name = run_params.iloc[0]['SN_name']
    # step = run_params.iloc[0]['n_steps']
    n_walkers = int(run_params.iloc[0]['n_walkers'])
    burn_in = int(run_params.iloc[0]['burn_in'])
    normalization = run_params.iloc[0]['normalization']
    LumThreshold = string_to_bool(run_params.iloc[0]['Tthreshold_lum'])
    data_dict = mcmc_snec.load_SN_data(data_type, SN_name)
    flat_sampler_path = os.path.join(result_path, 'flat_sampler.csv')
    sampler_df = pd.read_csv(flat_sampler_path,
                                        names=params,
                                        skiprows=(burn_in - 1) * (n_walkers))
    args = [SN_name, data_dict, sampler_df, ranges_dict, n_walkers, ax, normalization, LumThreshold, add_martinez]
    return args


def plot_result_fit(result_path, plot_types, ax, add_martinez=False):
    run_params_path = os.path.join(result_path, 'run_parameters.csv')
    run_params = pd.read_csv(run_params_path, index_col=0).T
    ranges_dict = import_ranges(run_params.columns.values, run_params)
    mcmc_snec.initialize_empty_models(ranges_dict)
    if 'lum' in plot_types:
        args = get_args_from_file(result_path, ax, 'lum', add_martinez)
        ax, log_likeli = plot_lum_with_fit(*args)
    if 'veloc' in plot_types:
        args = get_args_from_file(result_path, ax, 'veloc', add_martinez)
        ax, log_likeli =plot_veloc_with_fit(*args)
    if 'mag' in plot_types:
        args = get_args_from_file(result_path, ax, 'mag')
        ax, log_likeli =plot_mag_with_fit(*args)
    return ax, np.array(log_likeli)
