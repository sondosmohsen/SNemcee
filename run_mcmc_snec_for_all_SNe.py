'''

a
This script runs the mcmc_snec code based on the user-provided arguments:
Supernova (SN) name
number  of steps
number of walkers
output directory (usually would be the date)
by default, number of burn-in steps is 40% of total step number
'''

import numpy as np
import pandas as pd
import os
import mcmc_snec
import datetime
import sys, getopt

time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

'''
Lists specifying the range of parameters to allow the mcmc to explore,
based on the extent of the SNEC grid available
'''

Mzams_range = [9.0, 10.0, 11.0, 13.0, 15.0, 17.0]   # Progenitor ZAMS mass (solar mass)
Ni_range = [0.001, 0.02, 0.07, 0.12]# Ni mass (solar mass)
E_final_range = [0.1, 0.3, 0.5, 0.9, 1.3, 1.7]  # Explosion energy (foe=10^51 erg)
Mix_range = [2.0, 8.0]  # Ni mixing (solar mass)
R_range = [0, 500, 1000]
K_range = [0, 10, 30, 60]
T_range = [-5, 2]  # time shift from assumed day of explosion (days) T_range = [-5, 2]
parameter_ranges = {'Mzams': Mzams_range, 'Ni': Ni_range, 'E': E_final_range,
                    'R': R_range, 'K': K_range, 'Mix': Mix_range, 'S': [], 'T': T_range}
'''
The function loopy_snec_mcmc imports the SNEC_models for the SN, and then runs
mcmc_snec.py with the user-provided arguments: SN, step number, walker number 
as well as the run-type arguments determined by the main funciton:
*run_type: 'lum', 'veloc', 'mag', 'lum-veloc', 'mag-veloc', 'lum-mag', 'combined'
[which observations are used to calculate the fitting score in the likelihood function]
*csm: with, without, twostep, twostep_carryover
*LumThresh: True or False
*normalization: True or False
*nonuniform_priors: dictionary with keys being parameter names, and their values 
being more dictionaries with gaussian or polynomial distribution parameters. 
For example:
{'S': {'gaussian': {'mu': 1.0, 'sigma': sigma_S}}}
{'E': {'polynomial': np.poly1d(np.polyfit()) output}

The function outputs the MCMC SN_data, figures and running parameters used to
the output directory. Output 3.5es include:
1) Corner plot of posterior distributions for each parameter (all walkers, and steps after burn-in)        # TODO add this to the code
2) Fits of final models to the input observation SNEC_models for that SN (all walkers, and steps after burn-in)        # TODO add this to the code
3) Chain plots for each parameter        # TODO add this to the code
4) Flat_sampler.csv: MCMC sampler output table - a 2D table, where columns are the parameters, and rows are walkers and steps (by the hierarchy of: steps, and within them walkers)
5) final_results.csv: Summary statistics for the final posterior distributions for each parameter: mean, 16% percentile and 85% percentile.
6) run_parameters.csv: Record of the input arguments used in that run
'''

def make_res_dir(SN_name, run_type, csm, normalization, LumThreshold, output_dir, n_steps):
    filename = SN_name + '_' + run_type + '_csm-' + str(csm) \
               + '_normalized' + str(normalization) + '_TreshLum' + str(LumThreshold)
    res_dir = os.path.join('mcmc_results', output_dir, str(n_steps) + 'step', filename)
    if not os.path.exists(res_dir):
        if not os.path.exists(os.path.join('mcmc_results', output_dir, str(n_steps) + 'step')):
            if not os.path.exists(os.path.join('mcmc_results', output_dir)):
                os.mkdir(os.path.join('mcmc_results', output_dir))
            os.mkdir(os.path.join('mcmc_results', output_dir, str(n_steps) + 'step'))
        os.mkdir(res_dir)
    return res_dir


def loopy_snec_mcmc(SN_name, n_steps, burn_in, n_walkers, output_dir, parameter_ranges, run_type, csm, LumThreshold, normalization, nonuniform_priors):
    res_dir = make_res_dir(SN_name, run_type, csm, normalization, LumThreshold, output_dir, n_steps)
    S_range, sigma_S = mcmc_snec.S_prior(SN_name)
    SN_data_all = mcmc_snec.load_SN_data(run_type, SN_name)   #(sondos added this SN_data_all)
    parameter_ranges['S'] = S_range
    if nonuniform_priors is not None:
        nonuniform_priors = {'S': {'gaussian': {'mu': 1.0, 'sigma': sigma_S}}}
    '''
    # running code #
    '''
    if csm == 'with':
        sampler = mcmc_snec.emcee_fit_params(SN_name, res_dir, n_walkers, n_steps, burn_in,
                                                                      parameter_ranges, run_type, True, LumThreshold, normalization, nonuniform_priors)
        sampler_chain = sampler.chain
    elif csm == 'without' or csm == 'twostep' or csm == 'twostep-carryover':
        sampler = mcmc_snec.emcee_fit_params(SN_name, res_dir, n_walkers, n_steps, burn_in,
                                                                      parameter_ranges, run_type, False, LumThreshold, normalization, nonuniform_priors)
        sampler_chain = sampler.chain
        # The parameters R and K are not relevant for these cases, so they are filled with zeros by default.
        sampler_chain = np.insert(sampler_chain, 3, np.zeros((2, 1, 1)), axis=2)
    # make flat (2D) chain matrix
    s = list(sampler_chain.shape[1:])
    s[0] = np.prod(sampler_chain.shape[:2])
    sampler_chain_flat = sampler_chain.reshape(s)
    # make flat chain without the burn-in steps
    flat_sampler_no_burnin = sampler_chain_flat[n_walkers * burn_in:, :]
    final_step = n_steps - 1  # Python uses zero-based indexing, so the last step is n_steps - 1 (if it 500 steps, we will have: 0,1,.....499)


    if csm == 'twostep' or csm == 'twostep-carryover':
        # make a flat chain matrix as dataframe with column headers
        flat_sampler_no_burnin_df = pd.DataFrame(flat_sampler_no_burnin)
        flat_sampler_no_burnin_df.columns = ['Mzams', 'Ni', 'E', 'R', 'K', 'Mix', 'S', 'T']

        # get the positions of the walkers in the last step of the no-CSM stage
        last_walkers_noCSM = flat_sampler_no_burnin[-n_walkers:, :]
        # change the R and K of that last step to initial guesses of uniform distribution within their range
        last_walkers_noCSM[:, 3] = np.random.rand(n_walkers) * \
                                      (parameter_ranges['R'][-1] - parameter_ranges['R'][0]) + \
                                      parameter_ranges['R'][0]
        last_walkers_noCSM[:, 4] = np.random.rand(n_walkers) * \
                                      (parameter_ranges['K'][-1] - parameter_ranges['K'][0]) + \
                                      parameter_ranges['K'][0]
        # if we're doing carryover of the non-CSM results as priors to the second stage, extract their ditributions
        # in the last step and set as priors
        if csm == 'twostep':
            nonuniform_priors = None
        elif csm == 'twostep-carryover':
            nonuniform_priors = {}
            for param in ['Mzams', 'Ni', 'E', 'Mix', 'S', 'T']:
                nonuniform_priors[param] = {
                    'polynomial': mcmc_snec.polyfit_to_distribution(flat_sampler_no_burnin_df[param], res_dir)}

        # run the mcmc again, with walkers starting from their positions at the end of the first stage, and with or
        # without carryover as priors
        sampler_second_step = mcmc_snec.emcee_fit_params(SN_name, res_dir, n_walkers, n_steps, burn_in, parameter_ranges, #SN_data_all #(sondos added this SN_data_all)
                                                                                  run_type, True, LumThreshold, normalization, nonuniform_priors,
                                                                                  init_guesses=last_walkers_noCSM)
        # concatenate the second stage sampler to the first one
        #sampler_chain = np.concatenate((sampler_chain, sampler_second_step.chain), axis=1)
        #sampler_chain_flat = np.concatenate((sampler_chain_flat, sampler_second_step.get_chain(flat=True)), axis=0)
        #final_step = 2 * n_steps - 1
        sampler_chain = sampler_second_step.chain
        sampler_chain_flat = sampler_second_step.get_chain(flat=True)
        final_step = n_steps - 1

    #mcmc_snec.chain_plots(sampler_chain, parameter_ranges, res_dir, burn_in)
    #mcmc_snec.corner_plot(flat_sampler_no_burnin, parameter_ranges, res_dir)
    #mcmc_snec.plot_fit_with_data(sampler_chain,SN_data_all, parameter_ranges,
    #                             run_type, res_dir, SN_name, 0, LumThreshold, normalization)
    #mcmc_snec.plot_fit_with_data(sampler_chain,SN_data_all,parameter_ranges,
    #                             run_type, res_dir, SN_name, final_step, LumThreshold, normalization)
    mcmc_snec.save_param_results(sampler_chain, parameter_ranges, final_step, res_dir)
    np.savetxt(os.path.join(res_dir, 'flat_sampler.csv'), sampler_chain_flat, delimiter=",")

    print(sampler.chain.shape)

'''
Script can be run by calling run_mcmc_snec_for_all_SNe.py with user-provided
arguments: SN name, number of steps, number of steps and name for output directory) 
'''

def main(argv):
    SN_name = ''
    fitting_type = 'lum'
    n_steps = 500
    burn_in = 300
    n_walkers = 30
    output_dir = 'output_'+time_now
    csm = 'with'
    normalization = False
    LumThreshold = False
    nonuniform_priors = None
    arg_help = '{0} -S <SN name> -f <type of figure> -s <number of steps> -b <burn-in steps> -w <number of walkers> -o <output directory> -c <csm> -n <normalization> -l <luminosity threshold> -p <nonuniform priors>\n'\
               '\nargs:\n' \
               '-S SN name [required. for example: SN2017eaw]\n'\
               '-t fitting_type = lum, veloc, mag, lum-veloc, mag-veloc, lum-mag, lum-veloc-mag, combined [default: lum] \n' \
               '-s number of steps = <int> [default: 500]\n'\
               '-b burn-in steps = <int> [default: 300]\n'\
               '-w number of walkers = <int> [default: 30]\n' \
               '-o output directory name = <str> [default: output_<current time>]\n' \
               '-c csm = with, without, twostep, twostep-carryover [default: with] \n' \
               '-n normalization = True, False [default: False] \n' \
               '-l luminosity_threshold = True, False [default: False] \n'\
               '-p nonuniform_priors = None, <dictionary> [default: None] \n' \
               ''.format(argv[0])

    try:
        opts, args = getopt.getopt(argv[1:], "hS:t:s:b:w:o:c:n:l:p:", ["help", "SN=", 'fitting_type='
                                                                     "steps=", "burn_in=", "walkers=", "output_dir=",
                                                                     "csm=", "normalization=", "luminosity threshold=",
                                                                     "nonuniform_priors="])
    except:
        print(arg_help)
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-S", "--SN_name"):
            SN_name = arg
        elif opt in ("-t", "--fitting_type"):
            fitting_type = arg
        elif opt in ("-s", "--steps"):
            n_steps = int(arg)
        elif opt in ("-b", "--burn_in"):
            burn_in = int(arg)
        elif opt in ("-w", "--walkers"):
            n_walkers = int(arg)
        elif opt in ("-o", "--output_dir"):
            output_dir = arg
        elif opt in ("-c", "--csm"):
            csm = arg
        elif opt in ("-n", "--normalization"):
            normalization = arg
        elif opt in ("-l", "--luminosity_threshold"):
            LumThreshold = arg
        elif opt in ("-p", "--nonuniform_priors"):
            nonuniform_priors = arg

    print('SN name:', SN_name)
    print('fitting type:', fitting_type)
    print('steps:', n_steps)
    print('walkers:', n_walkers)
    print('burn in:', burn_in)
    print('output directory:', output_dir)
    print('csm:', csm)
    print('normalization:', normalization)
    print('luminosity threshold:', LumThreshold)
    print('nonuniform priors:', nonuniform_priors)

    loopy_snec_mcmc(SN_name, n_steps, burn_in, n_walkers, output_dir,
                    parameter_ranges, fitting_type, csm,
                    LumThreshold,
                    normalization,
                    nonuniform_priors)



if __name__ == "__main__":
    main(sys.argv)