# System
from copy import copy
from math import isnan, isinf
import re
# Dependencies from 3rd party
from molmass.molmass import Formula
import scipy.io
import numpy
import h5py
class inca_i():
    def __init__(self):
        self.fittedData=[];
        self.fittedFluxes=[];
        self.fittedFragments=[];
        self.fittedMeasuredFluxes=[];
        self.fittedMeasuredFragments=[];
        self.fittedMeasuredFluxResiduals=[];
        self.fittedMeasuredFragmentResiduals=[];
        self.simulationParameters=[];
    def import_isotopomerSimulationResults_INCA(self, simulation_id, filename, simulation_info, model_rxn_conversion_I=None):
        '''import results from a fluxomics simulation using INCA1.3
        INPUT:
        simulation_id = string, simulation_id
        filename = string, name of the matlab file
        simulation_info = {},
        model_rxn_conversion_I = optional {}, of INCA rxn_ids to model_rxn_ids (deprecated)
        Note: Please reference the model, fitdata, and simdata class structures in the INCA documentation
        for further information on the .mat file structure'''

        # extract information about the file
        import os, time
        from datetime import datetime
        from stat import ST_SIZE, ST_MTIME
        try:
            st = os.stat(filename)
        except IOError:
            print("failed to get information about", filename)
            return;
        else:
            file_size = st[ST_SIZE]
            simulation_dateAndTime_struct = time.localtime(st[ST_MTIME])
            simulation_dateAndTime = datetime.fromtimestamp(time.mktime(simulation_dateAndTime_struct))
        # determine if the simulation is a parallel labeling experiment, non-stationary, or both
        parallel = False;
        non_stationary = False;
        if len(simulation_info['experiment_id'])>1 or len(simulation_info['sample_name_abbreviation'])>1:
            parallel = True;
        if len(simulation_info['time_point'])>1:
            non_stationary = True;
        # extract out simulation data
        m = scipy.io.loadmat(filename)['m']; #model
        f = scipy.io.loadmat(filename)['f']; #fitdata
        s = scipy.io.loadmat(filename)['s']; #simdata
        # extract out model information (not currently recorded)
        m_ms_id = [];
        m_ms_on = [];
        m_ms_expt = [];
        for exp in m['expts']:
            exp_id = exp[0][0]['id'][0][0]
            for d in exp[0][0]['data_ms'][0]['id'][0]:
                m_ms_expt.append(exp_id);
                m_ms_id.append(d[0]);
                m_ms_on.append(bool(d[0][0]));
        # extract out simulation parameters (options)
        m_options = {
                    'cont_alpha':float(m['options'][0][0][0]['cont_alpha'][0][0][0]),
                    'cont_reltol':float(m['options'][0][0][0]['cont_reltol'][0][0][0]),
                    'cont_steps':float(m['options'][0][0][0]['cont_steps'][0][0][0]),
                    'fit_nudge':float(m['options'][0][0][0]['fit_nudge'][0][0][0]),
                    'fit_reinit':bool(m['options'][0][0][0]['fit_reinit'][0][0][0]),
                    'fit_reltol':float(m['options'][0][0][0]['fit_reltol'][0][0][0]),
                    'fit_starts':float(m['options'][0][0][0]['fit_starts'][0][0][0]),
                    'fit_tau':float(m['options'][0][0][0]['fit_tau'][0][0][0]),
                    'hpc_on':bool(m['options'][0][0][0]['hpc_on'][0][0][0]),
                    'int_maxstep':float(m['options'][0][0][0]['int_maxstep'][0][0][0]),
                    'int_reltol':float(m['options'][0][0][0]['int_reltol'][0][0][0]),
                    'int_senstol':float(m['options'][0][0][0]['int_senstol'][0][0][0]),
                    'int_timeout':float(m['options'][0][0][0]['int_timeout'][0][0][0]),
                    'int_tspan':float(m['options'][0][0][0]['int_tspan'][0][0][0]),
                    'ms_correct':bool(m['options'][0][0][0]['ms_correct'][0][0][0]),
                    'oed_crit':m['options'][0][0][0]['oed_crit'][0][0],
                    'oed_reinit':bool(m['options'][0][0][0]['oed_reinit'][0][0][0]),
                    'oed_tolf':float(m['options'][0][0][0]['oed_tolf'][0][0][0]),
                    'oed_tolx':float(m['options'][0][0][0]['oed_tolx'][0][0][0]),
                    'sim_more':bool(m['options'][0][0][0]['sim_more'][0][0][0]),
                    'sim_na':bool(m['options'][0][0][0]['sim_na'][0][0][0]),
                    'sim_sens':bool(m['options'][0][0][0]['sim_sens'][0][0][0]),
                    'sim_ss':bool(m['options'][0][0][0]['sim_ss'][0][0][0]),
                    'sim_tunit':m['options'][0][0][0]['sim_tunit'][0][0]
                    };
        try: m_options['hpc_mcr']=m['options'][0][0][0]['hpc_mcr'][0][0] #deprecated in INCA1.4
        except Exception as e:
            m_options['hpc_mcr']=float(m['options'][0][0][0]['hpc_bg'][0][0][0])
        try: m_options['hpc_serve']=m['options'][0][0][0]['hpc_serve'][0][0]#deprecated in INCA1.4
        except Exception as e:
            m_options['hpc_serve']=m['options'][0][0][0]['hpc_sched'][0][0]
        #if 'hpc_mcr' in m['options'][0][0][0]:
        #    m_options['hpc_mcr']=m['options'][0][0][0]['hpc_mcr'][0][0] #deprecated in INCA1.4
        #else:
        #    m_options['hpc_mcr']=float(m['options'][0][0][0]['hpc_bg'][0][0][0])
        #if 'hpc_serve' in m['options'][0][0][0]:
        #    m_options['hpc_serve']=m['options'][0][0][0]['hpc_serve'][0][0]#deprecated in INCA1.4
        #else:
        #    m_options['hpc_serve']=m['options'][0][0][0]['hpc_sched'][0][0]
        simulationParameters = [];
        m_options.update({'simulation_id':simulation_id,
		'simulation_dateAndTime':simulation_dateAndTime,
		'original_filename':filename,
                    'used_':True,
                    'comment_':None});
        simulationParameters.append(m_options);
        # extract out fit information
        fittedData = [];
        f_Echi2 = None;
        if not isnan(f['Echi2'][0][0][0][0]):
            if len(f['Echi2'][0][0][0])>1:
                f_Echi2 = [f['Echi2'][0][0][0][0],f['Echi2'][0][0][0][1]];
            else:
                f_Echi2 = [f['Echi2'][0][0][0][0]];
        f_alf = f['alf'][0][0][0][0];
        f_chi2 = f['chi2'][0][0][0][0];
        f_dof = int(f['dof'][0][0][0][0]);
        f_ = {'fitted_echi2':f_Echi2,
        'fitted_alf':f_alf,
        'fitted_chi2':f_chi2,
        'fitted_dof':f_dof};
        f_.update({'simulation_id':simulation_id,'simulation_dateAndTime':simulation_dateAndTime,
                    'used_':True,
                    'comment_':None})
        fittedData.append(f_);
        # extract out sum of the squared residuals of the fitted measurements
        f_mnt_id = [];
        f_mnt_sres = [];
        f_mnt_expt = [];
        f_mnt_type = []; #Flux or MS
        for d in f['mnt'][0][0][0]['id']:
            f_mnt_id.append(d[0]);
        for d in f['mnt'][0][0][0]['sres']:
            f_mnt_sres.append(float(d[0][0]));
        for d in f['mnt'][0][0][0]['expt']:
            f_mnt_expt.append(d[0]);
        for d in f['mnt'][0][0][0]['type']:
            f_mnt_type.append(d[0]);
        # seperate into appropriate table rows
        fittedMeasuredFluxes = [];
        fittedMeasuredFragments = [];
        for cnt,type in enumerate(f_mnt_type):
            if type=='Flux':
                if f_mnt_expt[cnt] in simulation_info['experiment_id']:
                    fittedMeasuredFluxes.append({'simulation_id':simulation_id,
                    'simulation_dateAndTime':simulation_dateAndTime,
                    'experiment_id':f_mnt_expt[cnt],
                    'sample_name_abbreviation':simulation_info['sample_name_abbreviation'][0],
                    'rxn_id':f_mnt_id[cnt],
                    'fitted_sres':f_mnt_sres[cnt],
                    'used_':True,
                    'comment_':None})
                elif f_mnt_expt[cnt] in simulation_info['sample_name_abbreviation']:
                    fittedMeasuredFluxes.append({'simulation_id':simulation_id,
                    'simulation_dateAndTime':simulation_dateAndTime,
                    'experiment_id':simulation_info['experiment_id'][0],
                    'sample_name_abbreviation':f_mnt_expt[cnt],
                    'rxn_id':f_mnt_id[cnt],
                    'fitted_sres':f_mnt_sres[cnt],
                    'used_':True,
                    'comment_':None})
            elif type=='MS':
                if f_mnt_expt[cnt] in simulation_info['experiment_id']:
                    fittedMeasuredFragments.append({'simulation_id':simulation_id,
                    'simulation_dateAndTime':simulation_dateAndTime,
                    'experiment_id':f_mnt_expt[cnt],
                    'sample_name_abbreviation':simulation_info['sample_name_abbreviation'][0],
                    'fragment_id':f_mnt_id[cnt],
                    'fitted_sres':f_mnt_sres[cnt],
                    'used_':True,
                    'comment_':None})
                elif f_mnt_expt[cnt] in simulation_info['sample_name_abbreviation']:
                    fittedMeasuredFragments.append({'simulation_id':simulation_id,
                    'simulation_dateAndTime':simulation_dateAndTime,
                    'experiment_id':simulation_info['experiment_id'][0],
                    'sample_name_abbreviation':f_mnt_expt[cnt],
                    'fragment_id':f_mnt_id[cnt],
                    'fitted_sres':f_mnt_sres[cnt],
                    'used_':True,
                    'comment_':None})
            else:
                print('type not recognized');
        # extract out the residuals of the fitted measurements
        f_mnt_res_val = [];
        f_mnt_res_fit = [];
        f_mnt_res_type = []; #Flux or MS
        f_mnt_res_id = [];
        f_mnt_res_std = [];
        f_mnt_res_time = [];
        f_mnt_res_expt = [];
        f_mnt_res_data = [];
        #f_mnt_res_esens = [];
        #f_mnt_res_msens = [];
        f_mnt_res_peak = [];
        for d in f['mnt'][0][0][0]['res']: 
            f_mnt_res_val.append(float(d[0][0]['val'][0][0]));
            f_mnt_res_fit.append(float(d[0][0]['fit'][0][0]));
            f_mnt_res_type.append(d[0][0]['type'][0]);
            f_mnt_res_id.append(d[0][0]['id'][0]);
            f_mnt_res_std.append(float(d[0][0]['std'][0][0]));
            #change default of time inf to 0
            if isinf(d[0][0]['time'][0][0]):
                f_mnt_res_time.append('0');
            else:
                f_mnt_res_time.append(str(d[0][0]['time'][0][0]));
            if d[0][0]['expt'][0] == 'Expt #1':
                f_mnt_res_expt.append(simulation_info['experiment_id'][0]);
            else:
                f_mnt_res_expt.append(d[0][0]['expt'][0]);
            f_mnt_res_data.append(d[0][0]['data'][0][0]);
            #f_mnt_res_esens.append(d[0][0]['esens'].data); #not needed, and matlab->python conversion has several bugs
            #f_mnt_res_msens.append(d[0][0]['msens'].data);
            if d[0][0]['peak']: f_mnt_res_peak.append(d[0][0]['peak'][0]);
            else: f_mnt_res_peak.append(None);
        # seperate into appropriate table rows
        fittedMeasuredFluxResiduals = [];
        fittedMeasuredFragmentResiduals = [];
        for cnt,type in enumerate(f_mnt_res_type):
            if type=='Flux':
                if f_mnt_res_expt[cnt] in simulation_info['experiment_id']:
                    fittedMeasuredFluxResiduals.append({'simulation_id':simulation_id,
                    'simulation_dateAndTime':simulation_dateAndTime,
                    'experiment_id':f_mnt_res_expt[cnt],
                    'sample_name_abbreviation':simulation_info['sample_name_abbreviation'][0],
                    'time_point':f_mnt_res_time[cnt],
                    'rxn_id':f_mnt_res_id[cnt],
                    'res_data':float(f_mnt_res_data[cnt]),
                    #'res_esens':f_mnt_res_esens[cnt],
                    'res_fit':float(f_mnt_res_fit[cnt]),
                    #'res_msens':f_mnt_res_msens[cnt],
                    'res_peak':f_mnt_res_peak[cnt],
                    'res_stdev':float(f_mnt_res_std[cnt]),
                    'res_val':float(f_mnt_res_val[cnt]),
                    'res_msens':None,
                    'res_esens':None,
                    'used_':True,
                    'comment_':None})
                elif f_mnt_res_expt[cnt] in simulation_info['sample_name_abbreviation']:
                    fittedMeasuredFluxResiduals.append({'simulation_id':simulation_id,
                    'simulation_dateAndTime':simulation_dateAndTime,
                    'experiment_id':simulation_info['experiment_id'][0],
                    'sample_name_abbreviation':f_mnt_res_expt[cnt],
                    'time_point':f_mnt_res_time[cnt],
                    'rxn_id':f_mnt_res_id[cnt],
                    'res_data':float(f_mnt_res_data[cnt]),
                    #'res_esens':f_mnt_res_esens[cnt],
                    'res_fit':float(f_mnt_res_fit[cnt]),
                    #'res_msens':f_mnt_res_msens[cnt],
                    'res_peak':f_mnt_res_peak[cnt],
                    'res_stdev':float(f_mnt_res_std[cnt]),
                    'res_val':float(f_mnt_res_val[cnt]),
                    'res_msens':None,
                    'res_esens':None,
                    'used_':True,
                    'comment_':None})
            elif type=='MS':
                # parse the id into fragment_id and mass
                fragment_string = f_mnt_res_id[cnt]
                fragment_string = re.sub('_DASH_','-',fragment_string)
                fragment_string = re.sub('_LPARANTHES_','[(]',fragment_string)
                fragment_string = re.sub('_RPARANTHES_','[)]',fragment_string)
                fragment_list = fragment_string.split('_');
                if not len(fragment_list)>5 or not ('MRM' in fragment_list or 'EPI' in fragment_list):
                    fragment_id = '_'.join([fragment_list[0],fragment_list[1],fragment_list[2]])
                    fragment_mass = Formula(fragment_list[2]).mass + float(fragment_list[3]);
                    time_point = fragment_list[4];
                else:
                    fragment_id = '_'.join([fragment_list[0],fragment_list[1],fragment_list[2],fragment_list[3]])
                    fragment_mass = Formula(fragment_list[2]).mass + float(fragment_list[4]);
                    time_point = fragment_list[5];
                #exp_id = fragment_list[5];
                #exp_id = fragment_list[6];
                fragment_id = re.sub('-','_DASH_',fragment_id)
                fragment_id = re.sub('[(]','_LPARANTHES_',fragment_id)
                fragment_id = re.sub('[)]','_RPARANTHES_',fragment_id)
                if f_mnt_res_expt[cnt] in simulation_info['experiment_id']:
                    fittedMeasuredFragmentResiduals.append({'simulation_id':simulation_id,
                    'simulation_dateAndTime':simulation_dateAndTime,
                    'experiment_id':f_mnt_res_expt[cnt],
                    'sample_name_abbreviation':simulation_info['sample_name_abbreviation'][0],
                    'time_point':f_mnt_res_time[cnt],
                    'fragment_id':fragment_id,
                    'fragment_mass':fragment_mass,
                    'res_data':float(f_mnt_res_data[cnt]),
                    #'res_esens':f_mnt_res_esens[cnt],
                    'res_fit':float(f_mnt_res_fit[cnt]),
                    #'res_msens':f_mnt_res_msens[cnt],
                    'res_peak':f_mnt_res_peak[cnt], #'res_peak':float(f_mnt_res_peak[cnt]),
                    'res_stdev':float(f_mnt_res_std[cnt]),
                    'res_val':float(f_mnt_res_val[cnt]),
                    'res_msens':None,
                    'res_esens':None,
                    'used_':True,
                    'comment_':None})
                elif f_mnt_res_expt[cnt] in simulation_info['sample_name_abbreviation']:
                    fittedMeasuredFragmentResiduals.append({'simulation_id':simulation_id,
                    'simulation_dateAndTime':simulation_dateAndTime,
                    'experiment_id':simulation_info['experiment_id'][0],
                    'sample_name_abbreviation':f_mnt_res_expt[cnt],
                    'time_point':f_mnt_res_time[cnt],
                    'fragment_id':fragment_id,
                    'fragment_mass':fragment_mass,
                    'res_data':float(f_mnt_res_data[cnt]),
                    #'res_esens':f_mnt_res_esens[cnt],
                    'res_fit':float(f_mnt_res_fit[cnt]),
                    #'res_msens':f_mnt_res_msens[cnt],
                    'res_peak':f_mnt_res_peak[cnt], #'res_peak':float(f_mnt_res_peak[cnt]),
                    'res_stdev':float(f_mnt_res_std[cnt]),
                    'res_val':float(f_mnt_res_val[cnt]),
                    'res_msens':None,
                    'res_esens':None,
                    'used_':True,
                    'comment_':None})
            else:
                print('type not recognized');
        # extract out the fitted parameters
        f_par_id = [];
        f_par_val = [];
        f_par_std = [];
        f_par_type = []; # 'Net flux' or 'Norm'
        f_par_lb = [];
        f_par_ub = [];
        f_par_unit = [];
        f_par_alf = [];
        f_par_chi2s = [];
        f_par_cor = [];
        f_par_cov = [];
        f_par_free = [];
        for d in f['par'][0][0][0]['id']:
            if 'Expt #1' in d[0]:
                id_str = d[0].astype('str')
                f_par_id.append(id_str.replace('Expt #1',simulation_info['experiment_id'][0]))
            else:
                f_par_id.append(d[0])
        # ensure that there are no negative values or infinite values
        for d in f['par'][0][0][0]['val']:
            if not d:
                f_par_val.append(0.0)
            #elif isnan(d[0][0]) or d[0][0]<1.0e-6:
            elif isnan(d[0][0]):
                f_par_val.append(0.0)
            #elif isinf(d[0][0]) or d[0][0]>1e3:
            elif isinf(d[0][0]):
                f_par_val.append(1.0e3)
            else:
                f_par_val.append(float(d[0][0]))
        for d in f['par'][0][0][0]['std']:
            if not d:
                f_par_std.append(0.0)
            elif isnan(d[0][0]):
                f_par_val.append(0.0)
            else:
                f_par_std.append(float(d[0][0]))
        for d in f['par'][0][0][0]['type']:
            f_par_type.append(d[0])
        #adjust the lb and ub to [0,1000];
        for cnt,d in enumerate(f['par'][0][0][0]['lb']):
            if not d:
                f_par_lb.append(0.0)
            #elif isnan(d[0][0]) or d[0][0]<1.0e-6:
            elif isnan(d[0][0]):
                f_par_lb.append(0.0)
            else:
                f_par_lb.append(float(d[0][0]))
        for d in f['par'][0][0][0]['ub']:
            if not d:
                f_par_ub.append(1.0e3)
            #elif isinf(d[0][0]) or isnan(d[0][0]) or d[0][0]>1.0e3:
            elif isinf(d[0][0]) or isnan(d[0][0]):
                f_par_ub.append(1.0e3)
            else:
                f_par_ub.append(float(d[0][0]))
        for d in f['par'][0][0][0]['unit']:
            if not d:
                #f_par_unit.append(None);
                #use default: mmol*gDCW-1*hr-1
                f_par_unit.append('mmol*gDCW-1*hr-1');
            else:
                f_par_unit.append(d[0])
        for d in f['par'][0][0][0]['alf']:
            f_par_alf.append(float(d[0][0]))
        for d in f['par'][0][0][0]['chi2s']:
            f_par_chi2s.append(d)
        for d in f['par'][0][0][0]['cor']:
            f_par_cor.append(d)
        for d in f['par'][0][0][0]['cov']:
            f_par_cov.append(d)
        for d in f['par'][0][0][0]['free']:
            f_par_free.append(bool(d[0][0]))
        # seperate into appropriate table rows
        fittedFluxes = [];
        fittedFragments = [];
        for cnt,type in enumerate(f_par_type):
            if type=='Net flux':
                fittedFluxes.append({'simulation_id':simulation_id,
                'simulation_dateAndTime':simulation_dateAndTime,
                'rxn_id':f_par_id[cnt],
                'flux':f_par_val[cnt],
                'flux_stdev':f_par_std[cnt],
                'flux_lb':f_par_lb[cnt],
                'flux_ub':f_par_ub[cnt],
                'flux_units':f_par_unit[cnt],
                'fit_alf':f_par_alf[cnt],
                'fit_chi2s':None,
                #'fit_chi2s':f_par_chi2s[cnt],
                'fit_cor':None,
                'fit_cov':None,
                #'fit_cor':f_par_cor[cnt],
                #'fit_cov':f_par_cov[cnt],
                'free':f_par_free[cnt],
                'used_':True,
                'comment_':None})
            elif type=='Norm':
                # parse the id 
                id_list = f_par_id[cnt].split(' ');
                expt = id_list[0];
                fragment_id = id_list[1];
                fragment_string = id_list[2];
                units = id_list[3];
                # parse the id into fragment_id and mass
                fragment_string = re.sub('_DASH_','-',fragment_string)
                fragment_string = re.sub('_LPARANTHES_','[(]',fragment_string)
                fragment_string = re.sub('_RPARANTHES_','[)]',fragment_string)
                fragment_list = fragment_string.split('_');
                if not len(fragment_list)>5 or not ('MRM' in fragment_list or 'EPI' in fragment_list):
                    fragment_mass = Formula(fragment_list[2]).mass + float(fragment_list[3]);
                    time_point = fragment_list[4];
                else:
                    fragment_mass = Formula(fragment_list[2]).mass + float(fragment_list[4]);
                    time_point = fragment_list[5];
                #fragment_id = '_'.join([fragment_list[0],fragment_list[1],fragment_list[2],fragment_list[3]])
                #fragment_id = re.sub('-','_DASH_',fragment_id)
                #fragment_id = re.sub('[(]','_LPARANTHES_',fragment_id)
                #fragment_id = re.sub('[)]','_RPARANTHES_',fragment_id)
                #expt_id = fragment_list[5];
                #expt_id = fragment_list[6];
                if expt in simulation_info['experiment_id']:
                    fittedFragments.append({'simulation_id':simulation_id,
                    'simulation_dateAndTime':simulation_dateAndTime,
                    'experiment_id':expt,
                    'sample_name_abbreviation':simulation_info['sample_name_abbreviation'][0],
                    'time_point':time_point,
                    'fragment_id':fragment_id,
                    'fragment_mass':fragment_mass,
                    'fit_val':f_par_val[cnt],
                    'fit_stdev':f_par_std[cnt],
                    'fit_units':units,
                    'fit_alf':f_par_alf[cnt],
                    'fit_cor':None,
                    'fit_cov':None,
                    #'fit_cor':f_par_cor[cnt],
                    #'fit_cov':f_par_cov[cnt],
                    'free':f_par_free[cnt],
                    'used_':True,
                    'comment_':None})
                elif expt in simulation_info['sample_name_abbreviation']:
                    fittedFragments.append({'simulation_id':simulation_id,
                    'simulation_dateAndTime':simulation_dateAndTime,
                    'experiment_id':simulation_info['experiment_id'][0],
                    'sample_name_abbreviation':expt,
                    'time_point':time_point,
                    'fragment_id':fragment_id,
                    'fragment_mass':fragment_mass,
                    'fit_val':f_par_val[cnt],
                    'fit_stdev':f_par_std[cnt],
                    'fit_units':units,
                    'fit_alf':f_par_alf[cnt],
                    'fit_cor':None,
                    'fit_cov':None,
                    #'fit_cor':f_par_cor[cnt],
                    #'fit_cov':f_par_cov[cnt],
                    'free':f_par_free[cnt],
                    'used_':True,
                    'comment_':None})
            else:
                print('type not recognized');
        # add data to the database
        self.fittedData=fittedData;
        self.fittedFluxes=fittedFluxes;
        self.fittedFragments=fittedFragments;
        self.fittedMeasuredFluxes=fittedMeasuredFluxes;
        self.fittedMeasuredFragments=fittedMeasuredFragments;
        self.fittedMeasuredFluxResiduals=fittedMeasuredFluxResiduals;
        self.fittedMeasuredFragmentResiduals=fittedMeasuredFragmentResiduals;
        self.simulationParameters=simulationParameters;