import os
from matplotlib import pyplot as plt
import plot_snec_fits as pltsn
import sys, getopt
import datetime
import numpy as np
import pandas as pd

time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")


def plot_single(fig_type, model_path, ax):
    ax, likeli = pltsn.plot_result_fit(model_path, fig_type, ax)
    ax.set_xlabel('Rest-Frame Days From Discovery', fontsize=14)
    if fig_type == 'lum':
        ax.set_ylabel('Log Bolometric Luminosity (erg/s)', fontsize=14)
    if fig_type == 'veloc':
        ax.set_ylabel('Expansion Velocity (km/s)', fontsize=14)
    if fig_type == 'mag':
        ax.set_ylabel('Absolute Magnitude', fontsize=14)
    plt.tight_layout()
    return ax, likeli


def composite_plot(SN_name, fig_type, fitting_type, csm, normalization, LumThreshold, results_dir, output_dir):
    fig_types = fig_type.split("-")
    model_name = SN_name + '_' + fitting_type + '_csm-' + csm + '_normalized' + str(normalization)\
                 + '_TreshLum' + str(LumThreshold)
    model_path = os.path.join(results_dir, model_name)
    num_subplots = len(fig_types)
    fig, axs = plt.subplots(num_subplots, figsize=(10, num_subplots*7))
    if num_subplots > 1:
        for i, f in enumerate(fig_types):
            axs[0], likeli = plot_single(f, model_path, axs[i])
            pltsn.add_likelihood_to_file(model_name, fig_type, likeli, output_dir)
    else:
        axs, likeli = plot_single(fig_type, model_path, axs)
        pltsn.add_likelihood_to_file(model_name, fig_type, likeli, output_dir)
    fig.savefig(os.path.join(output_dir, model_name + '_'+str(fig_type)+'_plot.png'))
    fig.savefig(os.path.join(output_dir, model_name + '_' + str(fig_type) + '_plot.pdf'))


def corner_plot(SN_name, fitting_type, csm, normalization, LumThreshold, results_dir, output_dir):
    model_name = SN_name + '_' + fitting_type + '_csm-' + csm + '_normalized' + str(normalization) \
                 + '_TreshLum' + str(LumThreshold)
    model_path = os.path.join(results_dir, model_name)
    pltsn.overlay_corner_plot([model_path], output_dir,
                                       [model_name], model_name)

def chain_plots(SN_name, fitting_type, csm, normalization, LumThreshold, results_dir, output_dir):
    model_name = SN_name + '_' + fitting_type + '_csm-' + csm + '_normalized' + str(normalization) \
                 + '_TreshLum' + str(LumThreshold)
    model_path = os.path.join(results_dir, model_name)
    pltsn.chain_plots(model_path, output_dir)


def lum_wCSM_vs_woCSM(SN_name, fitting_type, normalization, LumThreshold, results_dir, output_dir):
    fig, axs = plt.subplots(1, 2, sharey='row', figsize=(12, 6))    #default figsize=(20, 12))

    lum_csmTrue_name = SN_name + '_' + fitting_type + '_csm-with' + '_normalized' + str(normalization) \
                 + '_TreshLum' + str(LumThreshold)
    lum_csmTrue_path = os.path.join(results_dir, lum_csmTrue_name)

    lum_csmFalse_name = SN_name + '_' + fitting_type + '_csm-without' + '_normalized' + str(normalization) \
                       + '_TreshLum' + str(LumThreshold)
    lum_csmFalse_path = os.path.join(results_dir, lum_csmFalse_name)

    axs[0], likeli = pltsn.plot_result_fit(lum_csmTrue_path, 'lum', axs[0])
    pltsn.add_likelihood_to_file(lum_csmTrue_name, 'lum', likeli, output_dir)
    axs[1], likeli = pltsn.plot_result_fit(lum_csmFalse_path, 'lum', axs[1])
    pltsn.add_likelihood_to_file(lum_csmFalse_name, 'lum', likeli, output_dir)

    axs[0].set_ylabel('Log Bolometric Luminosity (erg/s)', fontsize=16)
    axs[0].set_xlabel('Rest-Frame Days From Discovery', fontsize=16)
    axs[1].set_xlabel('Rest-Frame Days From Discovery', fontsize=16)
    axs[0].yaxis.set_label_coords(-0.11, 0.5)
    axs[1].yaxis.set_label_coords(-0.11, 0.5)
    ################################################
    plt.ylim(41, 44)  # sondos added
    ################################################
    plt.tight_layout()
    #fig.savefig(os.path.join(output_dir, SN_name+'_lum_csm_comparison.png'))
    fig.savefig(os.path.join(output_dir, SN_name + '_lum_csm_comparison.pdf'))
    pltsn.overlay_corner_plot([lum_csmFalse_path, lum_csmTrue_path], output_dir,
                                       ['without CSM', 'with CSM'], SN_name + '_lum_csm_comparison')
    ################################################
    #axs[0].get_legend().remove()  # sondos added, the values should be in tables not  legend
    #axs[1].get_legend().remove()  # remove the values in the legend - the values should be in tables
    #axs[0].set_title('With CSM', fontsize=16)  # added title
    #axs[1].set_title('Without CSM', fontsize=16)  # added title
    ################################################
    #to see the values you need to delets this:
    # Add legends
    axs[0].scatter([], [], color='black', label="SN2017eaw")  # Add black dots to the legend
    axs[0].plot([], [], color='purple', label="fit with CSM")  # Add purple line to the legend
    axs[1].scatter([], [], color='black', label="SN2017eaw")  # Add black dots to the legend
    axs[1].plot([], [], color='purple', label="fit without CSM")  # Add purple line to the legend

    axs[0].legend(fontsize=18, loc='upper right')
    axs[1].legend(fontsize=18, loc='upper right')
    fig.savefig(os.path.join(output_dir, SN_name + '_lum_csm_comparison_final.pdf'))

    return fig


def lum_vs_lum_veloc_vs_lum_veloc_normalized(SN_name, csm, LumThreshold, results_dir, output_dir):
    fig, axs = plt.subplots(2, 3, sharey='row', figsize=(20, 10), gridspec_kw={'height_ratios': [2, 1.3]})

    lum_name = SN_name + '_lum_csm-' + csm + '_normalizedFalse_TreshLum' + str(LumThreshold)
    lum_path = os.path.join(results_dir, lum_name)
    lum_veloc_name = SN_name + '_lum-veloc_csm-' + csm + '_normalizedFalse_TreshLum' + str(LumThreshold)
    lum_veloc_path = os.path.join(results_dir, lum_veloc_name)
    lum_veloc_normalized_name = SN_name + '_lum-veloc_csm-' + csm + '_normalizedTrue_TreshLum' + str(LumThreshold)
    lum_veloc_normalized_path = os.path.join(results_dir, lum_veloc_normalized_name)

    axs[0, 0], likeli = pltsn.plot_result_fit(lum_path, 'lum', axs[0, 0])
    pltsn.add_likelihood_to_file(lum_name, 'lum', likeli, output_dir)
    axs[0, 1], likeli = pltsn.plot_result_fit(lum_veloc_path, 'lum', axs[0, 1])
    pltsn.add_likelihood_to_file(lum_veloc_name, 'lum', likeli, output_dir)
    axs[0, 2], likeli = pltsn.plot_result_fit(lum_veloc_normalized_path, 'lum', axs[0, 2])
    pltsn.add_likelihood_to_file(lum_veloc_normalized_name, 'lum', likeli, output_dir)
    axs[1, 0], likeli = pltsn.plot_result_fit(lum_path, 'veloc', axs[1, 0])
    pltsn.add_likelihood_to_file(lum_name, 'veloc', likeli, output_dir)
    axs[1, 1], likeli = pltsn.plot_result_fit(lum_veloc_path, 'veloc', axs[1, 1])
    pltsn.add_likelihood_to_file(lum_veloc_name, 'veloc', likeli, output_dir)
    axs[1, 2], likeli = pltsn.plot_result_fit(lum_veloc_normalized_path, 'veloc', axs[1, 2])
    pltsn.add_likelihood_to_file(lum_veloc_normalized_name, 'veloc', likeli, output_dir)

    axs[0, 0].set_ylabel('Log Bolometric Luminosity (erg/s)', fontsize=20)
    axs[1, 0].set_ylabel('Expansion Velocity (km/s)', fontsize=20)
    axs[1, 0].set_xlabel('Rest-Frame Days From Discovery', fontsize=20)
    axs[1, 1].set_xlabel('Rest-Frame Days From Discovery', fontsize=20)   #(lum_veloc)
    axs[1, 2].set_xlabel('Rest-Frame Days From Discovery', fontsize=20)    #(lum_veloc_normalized)
    axs[0, 0].yaxis.set_label_coords(-0.11, 0.5)
    axs[1, 0].yaxis.set_label_coords(-0.11, 0.5)

    axs[1, 0].set_ylim(3000, 11500), axs[1, 1].set_ylim(3000, 11500), axs[1, 2].set_ylim(3000, 11500)

    # Enlarge x and y ticks
    for ax_row in axs:
        for ax in ax_row:
            ax.tick_params(axis='both', which='major', labelsize=16)  # Enlarge tick labels for both axes
            ax.tick_params(axis='both', which='minor', labelsize=14)  # Optional for minor ticks

    plt.tight_layout()
    #fig.savefig(os.path.join(output_dir, SN_name + '_lum_veloc_comparison.png'))
    fig.savefig(os.path.join(output_dir, SN_name + '_lum_veloc_comparison.pdf'))

    pltsn.overlay_corner_plot([lum_path, lum_veloc_path, lum_veloc_normalized_path], output_dir,
                                       ['lum', 'lum+veloc, not normalized', 'lum+veloc, normalized'], SN_name + '_lum_veloc_comparison')
    ################################################
    axs[0,0].scatter([], [], color='black', label="SN2017eaw")
    axs[0,0].plot([], [], color='purple', label="lum")
    axs[0,1].scatter([], [], color='black', label="SN2017eaw")
    axs[0,1].plot([], [], color='purple', label="lum+veloc")
    axs[0,2].scatter([], [], color='black', label="SN2017eaw")
    axs[0,2].plot([], [], color='purple', label="lum+veloc, normalized")

    axs[0,0].legend(fontsize=18, loc='upper right')
    axs[0,1].legend(fontsize=18, loc='upper right')
    axs[0,2].legend(fontsize=18, loc='upper right')
    fig.savefig(os.path.join(output_dir, SN_name + '_lum_veloc_comparison_final.pdf'))

    return fig




def lum_veloc_vs_mag_veloc(SN_name, csm, normalization, LumThreshold, results_dir, output_dir):
    fig, axs = plt.subplots(3, 2, figsize=(14, 14), gridspec_kw={'height_ratios': [2,2,1.3]})

    lum_veloc_name = SN_name + '_lum-veloc_csm-' + csm + '_normalized' + str(normalization) \
                     + '_TreshLum' + str(LumThreshold)
    lum_veloc_path = os.path.join(results_dir, lum_veloc_name)
    mag_veloc_name = SN_name + '_mag-veloc_csm-' + csm + '_normalized' + str(normalization) \
                     + '_TreshLum' + str(LumThreshold)
    mag_veloc_path = os.path.join(results_dir, mag_veloc_name)

    axs[0, 0], likeli = pltsn.plot_result_fit(lum_veloc_path, 'lum', axs[0, 0])
    pltsn.add_likelihood_to_file(lum_veloc_name, 'lum', likeli, output_dir)
    axs[0, 1], likeli = pltsn.plot_result_fit(mag_veloc_path, 'lum', axs[0, 1])
    pltsn.add_likelihood_to_file(mag_veloc_name, 'lum', likeli, output_dir)
    axs[1, 0], likeli = pltsn.plot_result_fit(lum_veloc_path, 'mag', axs[1, 0])
    pltsn.add_likelihood_to_file(lum_veloc_name, 'mag', likeli, output_dir)
    axs[1, 1], likeli = pltsn.plot_result_fit(mag_veloc_path, 'mag', axs[1, 1])
    pltsn.add_likelihood_to_file(mag_veloc_name, 'mag', likeli, output_dir)
    axs[2, 0], likeli = pltsn.plot_result_fit(lum_veloc_path, 'veloc', axs[2, 0])
    pltsn.add_likelihood_to_file(lum_veloc_name, 'veloc', likeli, output_dir)
    axs[2, 1], likeli = pltsn.plot_result_fit(mag_veloc_path, 'veloc', axs[2, 1])
    pltsn.add_likelihood_to_file(mag_veloc_name, 'veloc', likeli, output_dir)

    axs[0, 0].set_ylabel('Log Bolometric Luminosity (erg/s)', fontsize=18)
    axs[1, 0].set_ylabel('Absolute Magnitude', fontsize=18)
    axs[2, 0].set_ylabel('Expansion Velocity (km/s)', fontsize=18)
    axs[2, 0].set_xlabel('Rest-Frame Days From Discovery', fontsize=18)
    axs[2, 1].set_xlabel('Rest-Frame Days From Discovery', fontsize=18)

    # Enlarge x and y ticks
    for ax_row in axs:
        for ax in ax_row:
            ax.tick_params(axis='both', which='major', labelsize=16)  # Enlarge tick labels for both axes
            ax.tick_params(axis='both', which='minor', labelsize=14)  # Optional for minor ticks

    plt.tight_layout()
    #fig.savefig(os.path.join(output_dir, SN_name + '_lum-veloc_mag-veloc_comparison_TreshLum' + str(LumThreshold) + '.png'))
    fig.savefig(
        os.path.join(output_dir, SN_name + '_lum-veloc_mag-veloc_comparison_TreshLum' + str(LumThreshold) + '.pdf'))

    pltsn.overlay_corner_plot([lum_veloc_path, mag_veloc_path], output_dir,
                                       ['lum+veloc', 'mag+veloc'],
                                       SN_name + '_lum-veloc_mag-veloc_comparison_TreshLum'+str(LumThreshold))
    ################################################
    axs[0, 0].scatter([], [], color='black', label="SN2017eaw")
    axs[0, 0].plot([], [], color='purple', label="lum+veloc")
    axs[0, 1].scatter([], [], color='black', label="SN2017eaw")
    axs[0, 1].plot([], [], color='purple', label="mag+veloc")

    axs[0, 0].legend(fontsize=18, loc='upper right')
    axs[0, 1].legend(fontsize=18, loc='upper right')
    fig.savefig(os.path.join(output_dir, SN_name + '_lum-veloc_mag-veloc_comparison_TreshLum' + str(LumThreshold) + 'final.pdf'))
    return fig




def lum_mag_veloc_Tthresh(SN_name, csm, normalization, results_dir, output_dir):
    fig, axs = plt.subplots(6, 2, figsize=(15, 21))

    lum_veloc_threshTrue_name = SN_name + '_lum-veloc_csm-' + csm + '_normalized' + str(normalization) \
                     + '_TreshLumTrue'
    lum_veloc_threshTrue_path = os.path.join(results_dir, lum_veloc_threshTrue_name)
    mag_veloc_threshTrue_name = SN_name + '_mag-veloc_csm-' + csm + '_normalized' + str(normalization) \
                     + '_TreshLumTrue'
    mag_veloc_threshTrue_path = os.path.join(results_dir, mag_veloc_threshTrue_name)

    lum_veloc_threshFalse_name = SN_name + '_lum-veloc_csm-' + csm + '_normalized' + str(normalization) \
                                + '_TreshLumFalse'
    lum_veloc_threshFalse_path = os.path.join(results_dir, lum_veloc_threshFalse_name)
    mag_veloc_threshFalse_name = SN_name + '_mag-veloc_csm-' + csm + '_normalized' + str(normalization) \
                                + '_TreshLumFalse'
    mag_veloc_threshFalse_path = os.path.join(results_dir, mag_veloc_threshFalse_name)


    axs[0, 0], likeli = pltsn.plot_result_fit(lum_veloc_threshTrue_path, 'lum', axs[0, 0])
    pltsn.add_likelihood_to_file(lum_veloc_threshTrue_name, 'lum', likeli, output_dir)
    axs[0, 1], likeli = pltsn.plot_result_fit(mag_veloc_threshTrue_path, 'lum', axs[0, 1])
    pltsn.add_likelihood_to_file(mag_veloc_threshTrue_name, 'lum', likeli, output_dir)
    axs[1, 0], likeli = pltsn.plot_result_fit(lum_veloc_threshTrue_path, 'mag', axs[1, 0])
    pltsn.add_likelihood_to_file(lum_veloc_threshTrue_name, 'mag', likeli, output_dir)
    axs[1, 1], likeli = pltsn.plot_result_fit(mag_veloc_threshTrue_path, 'mag', axs[1, 1])
    pltsn.add_likelihood_to_file(mag_veloc_threshTrue_name, 'mag', likeli, output_dir)
    axs[2, 0], likeli = pltsn.plot_result_fit(lum_veloc_threshTrue_path, 'veloc', axs[2, 0])
    pltsn.add_likelihood_to_file(lum_veloc_threshTrue_name, 'veloc', likeli, output_dir)
    axs[2, 1], likeli = pltsn.plot_result_fit(mag_veloc_threshTrue_path, 'veloc', axs[2, 1])
    pltsn.add_likelihood_to_file(mag_veloc_threshTrue_name, 'veloc', likeli, output_dir)

    axs[3, 0], likeli = pltsn.plot_result_fit(lum_veloc_threshFalse_path, 'lum', axs[3, 0])
    pltsn.add_likelihood_to_file(lum_veloc_threshFalse_name, 'lum', likeli, output_dir)
    axs[3, 1], likeli = pltsn.plot_result_fit(mag_veloc_threshFalse_path, 'lum', axs[3, 1])
    pltsn.add_likelihood_to_file(mag_veloc_threshFalse_name, 'lum', likeli, output_dir)
    axs[4, 0], likeli = pltsn.plot_result_fit(lum_veloc_threshFalse_path, 'mag', axs[4, 0])
    pltsn.add_likelihood_to_file(lum_veloc_threshFalse_name, 'mag', likeli, output_dir)
    axs[4, 1], likeli = pltsn.plot_result_fit(mag_veloc_threshFalse_path, 'mag', axs[4, 1])
    pltsn.add_likelihood_to_file(mag_veloc_threshFalse_name, 'mag', likeli, output_dir)
    axs[5, 0], likeli = pltsn.plot_result_fit(lum_veloc_threshFalse_path, 'veloc', axs[5, 0])
    pltsn.add_likelihood_to_file(lum_veloc_threshFalse_name, 'veloc', likeli, output_dir)
    axs[5, 1], likeli = pltsn.plot_result_fit(mag_veloc_threshFalse_path, 'veloc', axs[5, 1])
    pltsn.add_likelihood_to_file(mag_veloc_threshFalse_name, 'veloc', likeli, output_dir)

    axs[0, 0].set_ylabel('Log Bolometric Luminosity (erg/s)', fontsize=12)
    axs[1, 0].set_ylabel('Absolute Magnitude', fontsize=12)
    axs[2, 0].set_ylabel('Expansion Velocity (km/s)', fontsize=12)
    axs[2, 0].set_xlabel('Rest-Frame Days From Discovery', fontsize=12)
    axs[2, 1].set_xlabel('Rest-Frame Days From Discovery', fontsize=12)

    axs[3, 0].set_ylabel('Log Bolometric Luminosity (erg/s)', fontsize=12)
    axs[4, 0].set_ylabel('Absolute Magnitude', fontsize=12)
    axs[5, 0].set_ylabel('Expansion Velocity (km/s)', fontsize=12)
    axs[5, 0].set_xlabel('Rest-Frame Days From Discovery', fontsize=12)
    axs[5, 1].set_xlabel('Rest-Frame Days From Discovery', fontsize=12)


    plt.tight_layout()
    #fig.savefig(os.path.join(output_dir, SN_name + '_lum-veloc_mag-veloc_comparison_TreshLumComparison.png'))
    fig.savefig(
        os.path.join(output_dir, SN_name + '_lum-veloc_mag-veloc_comparison_TreshLumComparison.pdf'))

    return fig



def lum_veloc_onestep_vs_twostep(SN_name, normalization, LumThreshold, results_dir, output_dir):
    fig, axs = plt.subplots(2, 3, sharey='row', figsize=(20, 10), gridspec_kw={'height_ratios': [2, 1.3]})
    onestep_name = SN_name + '_lum-veloc_csm-with_normalized' \
                   + str(normalization) + '_TreshLum' + str(LumThreshold)
    onestep_path = os.path.join(results_dir, onestep_name)
    twostep_priorNone_name = SN_name + '_lum-veloc_csm-twostep_normalized' \
                             + str(normalization) + '_TreshLum' + str(LumThreshold)
    twostep_priorNone_path = os.path.join(results_dir, twostep_priorNone_name)
    twostep_priorTrue_name = SN_name + '_lum-veloc_csm-twostep-carryover_normalized' \
                             + str(normalization) + '_TreshLum' + str(LumThreshold)
    twostep_priorTrue_path = os.path.join(results_dir, twostep_priorTrue_name)

    axs[0, 0], likeli = pltsn.plot_result_fit(onestep_path, 'lum', axs[0, 0])
    pltsn.add_likelihood_to_file(onestep_name, 'lum', likeli, output_dir)
    axs[0, 1], likeli = pltsn.plot_result_fit(twostep_priorNone_path, 'lum', axs[0, 1])
    pltsn.add_likelihood_to_file(twostep_priorNone_name, 'lum', likeli, output_dir)
    axs[0, 2], likeli = pltsn.plot_result_fit(twostep_priorTrue_path, 'lum', axs[0, 2])
    pltsn.add_likelihood_to_file(twostep_priorTrue_name, 'lum', likeli, output_dir)
    axs[1, 0], likeli = pltsn.plot_result_fit(onestep_path, 'veloc', axs[1, 0])
    pltsn.add_likelihood_to_file(onestep_name, 'veloc', likeli, output_dir)
    axs[1, 1], likeli = pltsn.plot_result_fit(twostep_priorNone_path, 'veloc', axs[1, 1])
    pltsn.add_likelihood_to_file(twostep_priorNone_name, 'veloc', likeli, output_dir)
    axs[1, 2], likeli = pltsn.plot_result_fit(twostep_priorTrue_path, 'veloc', axs[1, 2])
    pltsn.add_likelihood_to_file(twostep_priorTrue_name, 'veloc', likeli, output_dir)

    axs[1, 0].set_ylim(3000, 11500), axs[1, 1].set_ylim(3000, 11500), axs[1, 2].set_ylim(3000, 11500)

    axs[0, 0].set_ylabel('Log Bolometric Luminosity (erg/s)', fontsize=20)
    axs[1, 0].set_ylabel('Expansion Velocity (km/s)', fontsize=20)
    axs[1, 0].set_xlabel('Rest-Frame Days From Discovery', fontsize=20)
    axs[1, 1].set_xlabel('Rest-Frame Days From Discovery', fontsize=20)
    axs[1, 2].set_xlabel('Rest-Frame Days From Discovery', fontsize=20)
    axs[0, 0].yaxis.set_label_coords(-0.11, 0.5)
    axs[1, 0].yaxis.set_label_coords(-0.11, 0.5)

    # Enlarge x and y ticks
    for ax_row in axs:
        for ax in ax_row:
            ax.tick_params(axis='both', which='major', labelsize=16)  # Enlarge tick labels for both axes
            ax.tick_params(axis='both', which='minor', labelsize=14)  # Optional for minor ticks

    plt.tight_layout()
    #fig.savefig(os.path.join(output_dir, SN_name + '_onestep_twostep_comparison.png'))
    fig.savefig(os.path.join(output_dir, SN_name + '_onestep_twostep_comparison.pdf'))

    pltsn.overlay_corner_plot([onestep_path, twostep_priorNone_path, twostep_priorTrue_path],
                                       output_dir,
                                       ['one step', 'two step, no priors carryover', 'two step, with priors carryover'],
                                       SN_name + '_onestep_twostep_comparison')
    ################################################
    axs[0, 0].scatter([], [], color='black', label="SN2017eaw")
    axs[0, 0].plot([], [], color='purple', label="One step")
    axs[0, 1].scatter([], [], color='black', label="SN2017eaw")
    axs[0, 1].plot([], [], color='purple', label="Two steps ")
    axs[0, 2].scatter([], [], color='black', label="SN2017eaw")
    axs[0, 2].plot([], [], color='purple', label="Two steps with priors carryover")

    axs[0, 0].legend(fontsize=18, loc='upper right')
    axs[0, 1].legend(fontsize=18, loc='upper right')
    axs[0, 2].legend(fontsize=18, loc='upper right')
    fig.savefig(os.path.join(output_dir, SN_name + '_onestep_twostep_comparison_final.pdf'))
    return fig


def lum_veloc_vs_martinez1(SN_name, results_dir, output_dir):   #(old version (just with csm))
    fig, axs = plt.subplots(2, 1, sharey='row', figsize=(8, 12))   #default figsize=(20, 16)
    lum_veloc_name = SN_name + '_lum-veloc_csm-with_normalizedFalse_TreshLumFalse'
    lum_veloc_path = os.path.join(results_dir, lum_veloc_name)
    pltsn.plot_result_fit(lum_veloc_path, 'lum', axs[0], add_martinez=True)
    pltsn.plot_result_fit(lum_veloc_path, 'veloc', axs[1], add_martinez=True)
    axs[0].set_ylabel('Log Bolometric Luminosity (erg/s)', fontsize=22) #default fontsize=14
    axs[1].set_ylabel('Expansion Velocity (km/s)', fontsize=22)
    axs[1].set_xlabel('Rest-Frame Days From Discovery', fontsize=22)
    axs[0, 0].yaxis.set_label_coords(-0.11, 0.5)  # Adjust x and y positions
    axs[1, 0].yaxis.set_label_coords(-0.11, 0.5)  # Same x position for consistency
    # Enlarge x and y ticks
    for ax in axs:
        ax.tick_params(axis='both', which='major', labelsize=16)  # Enlarge tick labels for both axes
        ax.tick_params(axis='both', which='minor', labelsize=12)  # Optional for minor ticks
    plt.tight_layout()
    fig.savefig(os.path.join('figures',output_dir, SN_name + '_martinez_comparison.pdf'))   #sondos added: ('figures', to save the figure into the same direc
    ################################################
    axs[0].scatter([], [], color='black', label=SN_name)
    axs[0].plot([], [], color='purple', label="SNEmcee")
    axs[0].plot([], [], color='orange', label="Martinez values on SNEmcee ")
    axs[0].legend(fontsize=18, loc='upper right')
    axs[1].get_legend().remove()
    fig.savefig(os.path.join('figures',output_dir, SN_name + '_martinez_comparison_final.pdf'))
    ################################################ adding bestfit models of martinez
    # Adding Martinez best-fit models
    models_path = os.path.join('results', 'martinez_bestfit_models', f"{SN_name}.csv")
    models = pd.read_csv(models_path)

    if 'VPH' in models.columns:  # Convert velocity from cm/s to km/s
        models['VPH'] = models['VPH'] / 1e5  # Convert from cm/s to km/s

    axs[0].plot(models['TIME'], models['LOGL'], color='green', label='Martinez best-fit models')

    velocity_data = models[models['TIME'] <= 125]
    axs[1].plot(velocity_data['TIME'], velocity_data['VPH'], color='green', label='Martinez Best-fit')

    # Update the legends
    axs[0].legend(fontsize=18, loc='upper right')

    fig.savefig(os.path.join('figures', output_dir, SN_name + '_martinez_comparison_with_bestfit.pdf'))

    return fig


def lum_veloc_vs_martinez2(SN_name, results_dir, output_dir):
    fig, axs = plt.subplots(2, 2, sharey='row', figsize=(12, 8), gridspec_kw={'height_ratios': [2, 1.5]})

    # File paths for with-CSM and without-CSM
    lum_veloc_with_csm_name = SN_name + '_lum-veloc_csm-with_normalizedFalse_TreshLumFalse'
    lum_veloc_with_csm_path = os.path.join(results_dir, lum_veloc_with_csm_name)

    lum_veloc_without_csm_name = SN_name + '_lum-veloc_csm-without_normalizedFalse_TreshLumFalse'
    lum_veloc_without_csm_path = os.path.join(results_dir, lum_veloc_without_csm_name)

    # Plotting bolometric luminosity
    pltsn.plot_result_fit(lum_veloc_with_csm_path, 'lum', axs[0, 0], add_martinez=True)
    pltsn.plot_result_fit(lum_veloc_without_csm_path, 'lum', axs[0, 1], add_martinez=True)

    # Plotting velocities
    pltsn.plot_result_fit(lum_veloc_with_csm_path, 'veloc', axs[1, 0], add_martinez=True)
    pltsn.plot_result_fit(lum_veloc_without_csm_path, 'veloc', axs[1, 1], add_martinez=True)

    axs[1, 0].set_ylim(3000, 11000)
    axs[1, 1].set_ylim(3000, 11000)

    # Labels for each subplot
    axs[0, 0].set_ylabel('Log Bolometric Luminosity (erg/s)', fontsize=16)
    axs[1, 0].set_ylabel('Expansion Velocity (km/s)', fontsize=16)
    axs[1, 0].set_xlabel('Rest-Frame Days From Discovery', fontsize=16)
    axs[1, 1].set_xlabel('Rest-Frame Days From Discovery', fontsize=16)

    axs[0, 0].yaxis.set_label_coords(-0.11, 0.5)  # Adjust x and y positions
    axs[1, 0].yaxis.set_label_coords(-0.11, 0.5)  # Same x position for consistency

    # Enlarge x and y ticks for all subplots
    for ax_row in axs:
        for ax in ax_row:
            ax.tick_params(axis='both', which='major', labelsize=12)
            ax.tick_params(axis='both', which='minor', labelsize=12)

    # Adding legends
    axs[0, 0].scatter([], [], color='black', label=SN_name)
    axs[0, 0].plot([], [], color='purple', label="SNEmcee with CSM")
    axs[0, 0].plot([], [], color='orange', label="Martinez values on SNEmcee")
    axs[0, 0].legend(fontsize=14, loc='upper right')
    axs[0, 1].scatter([], [], color='black', label=SN_name)
    axs[0, 1].plot([], [], color='purple', label="SNEmcee without CSM")
    axs[0, 1].plot([], [], color='orange', label="Martinez values on SNEmcee")
    axs[0, 1].legend(fontsize=14, loc='upper right')

    # Adding Martinez best-fit models
    models_path = os.path.join('results', 'martinez_bestfit_models', f"{SN_name}.csv")
    models = pd.read_csv(models_path)

    if 'VPH' in models.columns:
        models['VPH'] = models['VPH'] / 1e5  # Convert from cm/s to km/s

    axs[0,0].plot(models['TIME'], models['LOGL'], color='green', label='Martinez best-fit models')
    axs[0,1].plot(models['TIME'], models['LOGL'], color='green', label='Martinez best-fit models')
    axs[0, 0].legend(fontsize=14, loc='upper right')
    velocity_data = models[models['TIME'] <= 125]
    axs[1, 0].plot(velocity_data['TIME'], velocity_data['VPH'], color='green', label='Martinez Best-fit')
    axs[1, 1].plot(velocity_data['TIME'], velocity_data['VPH'], color='green', label='Martinez Best-fit')
    axs[0, 1].legend(fontsize=14, loc='upper right')
    plt.tight_layout()
    fig.savefig(os.path.join('figures', output_dir, SN_name + '_martinez_comparison_with_bestfit1.pdf'))
    axs[1, 0].get_legend().remove()
    axs[1, 1].get_legend().remove()
    fig.savefig(os.path.join('figures', output_dir, SN_name + '_martinez_comparison_with_bestfit.pdf'))

    return fig


def lum_veloc_vs_martinez(SN_name, results_dir, output_dir):
    import os
    import pandas as pd
    import matplotlib.pyplot as plt

    # Create subplots: 2 stacked vertically
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(8, 12), gridspec_kw={'height_ratios': [2, 1.5]})

    # File paths for with-CSM and without-CSM
    lum_veloc_with_csm_name = SN_name + '_lum-veloc_csm-with_normalizedFalse_TreshLumFalse'
    lum_veloc_with_csm_path = os.path.join(results_dir, lum_veloc_with_csm_name)

    lum_veloc_without_csm_name = SN_name + '_lum-veloc_csm-without_normalizedFalse_TreshLumFalse'
    lum_veloc_without_csm_path = os.path.join(results_dir, lum_veloc_without_csm_name)

    # Extract data from plot_result_fit (assuming it returns the line object or data)
    pltsn.plot_result_fit(lum_veloc_with_csm_path, 'lum', axs[0], add_martinez=True)
    pltsn.plot_result_fit(lum_veloc_without_csm_path, 'lum', axs[0], add_martinez=True)
    pltsn.plot_result_fit(lum_veloc_with_csm_path, 'veloc', axs[1], add_martinez=True)
    pltsn.plot_result_fit(lum_veloc_without_csm_path, 'veloc', axs[1], add_martinez=True)

    # Assign colors to ensure correct identification
    with_csm_lum_line = axs[0].lines[-4]  # Assuming the order of addition is predictable
    without_csm_lum_line = axs[0].lines[-3]
    with_csm_vel_line = axs[1].lines[-4]
    without_csm_vel_line = axs[1].lines[-3]

    # Set colors for with-CSM and without-CSM
    with_csm_lum_line.set_color('red')
    without_csm_lum_line.set_color('blue')
    with_csm_vel_line.set_color('red')
    without_csm_vel_line.set_color('blue')

    # Adjust velocity limits for consistency
    axs[1].set_ylim(3000, 11000)

    # Add labels
    axs[0].set_ylabel('Log Bolometric Luminosity (erg/s)', fontsize=20)
    axs[1].set_ylabel('Expansion Velocity (km/s)', fontsize=20)
    axs[1].set_xlabel('Rest-Frame Days From Discovery', fontsize=20)

    axs[0].yaxis.set_label_coords(-0.11, 0.5)  # Adjust x and y positions for top plot
    axs[1].yaxis.set_label_coords(-0.11, 0.5)  # Same x position for bottom plot

    # Enlarge x and y ticks for both subplots
    for ax in axs:
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.tick_params(axis='both', which='minor', labelsize=14)

    # Add legends
    axs[0].scatter([], [], color='black', label=SN_name)
    axs[0].plot([], [], color='red', label="SNEmcee with CSM")
    axs[0].plot([], [], color='blue', label="SNEmcee without CSM")
    axs[0].plot([], [], color='orange', label="Martinez values on SNEmcee")
    axs[0].legend(fontsize=14, loc='upper right')

    # Adding Martinez best-fit models
    models_path = os.path.join('results', 'martinez_bestfit_models', f"{SN_name}.csv")
    models = pd.read_csv(models_path)

    if 'VPH' in models.columns:  # Convert velocity from cm/s to km/s
        models['VPH'] = models['VPH'] / 1e5  # Convert to km/s

    axs[0].plot(models['TIME'], models['LOGL'], color='green', label='Martinez Best-fit Models')
    velocity_data = models[models['TIME'] <= 125]
    axs[1].plot(velocity_data['TIME'], velocity_data['VPH'], color='green', label='Martinez Best-fit')

    axs[0].legend(fontsize=14, loc='upper right')

    # Save intermediate figure
    fig.savefig(os.path.join('figures', output_dir, SN_name + '_martinez_comparison.pdf'))

    # Add final touches and save the final figure
    axs[1].get_legend().remove()  # Remove bottom legend
    fig.savefig(os.path.join('figures', output_dir, SN_name + '_martinez_comparison_with_bestfit.pdf'))

    return fig







def main(argv):
    SN_name = ''
    n_steps = 500
    type_fig = 'lum'
    output_dir = 'output_' + time_now
    fitting_type = 'lum'
    csm = 'without'
    normalization = False
    LumThreshold = False

    arg_help = '{0} -S <SN name> -s <number of steps> -f <figure_type> -o <output directory> -t <fitting type> -c <csm> -n <normalization> -l <luminosity threshold>\n'\
               '\nargs:\n' \
               '-S SN name [required. for example: SN2017eaw]\n' \
               '-s number of steps = <int> [default: 500]\n' \
               '-f figure type = single fit-type plots: lum/veloc/mag [or any combination of those, separated by -] ' \
               'OR walkers OR corner OR comparison plots: ' \
               'csm_comparison/lum-mag-veloc_comparison/lum-veloc-normalized_comparison/lum-veloc-twostep_comparison/' \
               'martinez_comparison [required]\n' \
               '-o output directory name = <str> [default: output_<current time>]\n' \
               '-t fitting_type = lum, veloc, mag, lum-veloc, mag-veloc, lum-mag, lum-veloc-mag [default: False] \n' \
               '-c csm = with, without, twostep, twostep-carryover [default: False] \n' \
               '-n normalization = True, False [default: False] \n' \
               '-l luminosity_threshold = True, False [default: False] \n' \
               ''.format(argv[0])
    try:
        opts, args = getopt.getopt(argv[1:], "hS:s:f:o:t:c:l:n:", ["help", "SN=", "steps=", "type of figure=", "output_dir=",
                                                                   "fitting type=", "csm=", "normalization=", "luminosity threshold="])
    except:
        print(arg_help)
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-S", "--SN_name"):
            SN_name = arg
        elif opt in ("-s", "--steps"):
            n_steps = int(arg)
        elif opt in ("-f", "--figure_type"):
            type_fig = arg
        elif opt in ("-o", "--output_dir"):
            output_dir = arg
        elif opt in ("-t", "--fitting_type"):
            fitting_type = arg
        elif opt in ("-c", "--csm"):
            csm = arg
        elif opt in ("-n", "--normalization"):
            normalization = arg
        elif opt in ("-l", "--luminosity_threshold"):
            LumThreshold = arg
    print('SN name:', SN_name)
    print('steps:', n_steps)
    print('output directory:', output_dir)
    print('type of figure:', type_fig)
    print('fitting type:', fitting_type)
    print('csm:', csm)
    print('normalization:', normalization)
    print('luminosity threshold:', LumThreshold)

    res_dir = os.path.join('mcmc_results', output_dir, str(n_steps) + 'step')
    step_dir = os.path.join('figures', output_dir, str(n_steps) + 'step')
    if not os.path.exists(step_dir):
        if not os.path.exists(os.path.join('figures', output_dir)):
            os.mkdir(os.path.join('figures', output_dir))
        os.mkdir(step_dir)

    if 'comparison' not in type_fig:
        if type_fig == 'corner':
            corner_plot(SN_name, fitting_type, csm, normalization, LumThreshold, res_dir, step_dir)
        elif type_fig == 'walkers':
            chain_plots(SN_name, fitting_type, csm, normalization, LumThreshold, res_dir, step_dir)
        else:
            composite_plot(SN_name, type_fig, fitting_type, csm, normalization, LumThreshold, res_dir, step_dir)
    elif type_fig == 'csm_comparison':
        lum_wCSM_vs_woCSM(SN_name, fitting_type, normalization, LumThreshold, res_dir, step_dir)
    elif type_fig == 'lum-mag-veloc_comparison':
        lum_veloc_vs_mag_veloc(SN_name, csm, normalization, LumThreshold, res_dir, step_dir)
    elif type_fig == 'lum-mag-veloc-Tthresh_comparison':       #sondos changed Tthersh to Tthresh
        lum_mag_veloc_Tthresh(SN_name, csm, normalization, res_dir, step_dir)  # sondos chaged to lum_mag_veloc_Tthresh
        #lum_veloc_vs_mag_veloc(SN_name, csm, normalization, res_dir, step_dir)      #sondos added LumThreshold
    elif type_fig == 'lum-veloc-normalized_comparison':
        lum_vs_lum_veloc_vs_lum_veloc_normalized(SN_name, csm, LumThreshold, res_dir, step_dir)
    elif type_fig == 'lum-veloc-twostep_comparison':
        lum_veloc_onestep_vs_twostep(SN_name, normalization, LumThreshold, res_dir, step_dir)
    elif type_fig == 'martinez_comparison':
        lum_veloc_vs_martinez(SN_name, res_dir, output_dir)
    else:
        print('figure type (-f) should be csm_comparison, lum-mag-veloc_comparison, lum-veloc-normalized_comparison, '
              'lum-veloc-twostep_comparison,lum-mag-veloc-Tthresh_comparison, martinez_comparison, or any combination of lum, mag and veloc separated by'
              ' a dash (-)')

if __name__ == "__main__":
    main(sys.argv)
