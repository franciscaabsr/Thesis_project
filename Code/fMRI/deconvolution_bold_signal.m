%% Script for computing the hemodynamic response delay for each ROI

% considering rsHRF_estimation_temporal_basis.m so as to use a standard
% hemodynamic response, bf, including temporal and dispersion derivatives

% can also be done with smoothed FIR, making fit to the data so as to
% obtain the impulse response

% Code adapted from demo examples of rsHRF: a Toolbox for Resting State HRF
% Deconvolution and Connectivity Analysis by Wu et al.

% demo code for voxel-wise HRF deconvolution
% From NIFTI image (resting state fMRI data) to NIFTI image (HRF parameters).
% Guo-Rong Wu, gronwu@gmail.com, SWU, 2019.10.30
% Reference: Wu, G.; Liao, W.; Stramaglia, S.; Ding, J.; Chen, H. & Marinazzo, D..
% A blind deconvolution approach to recover effective connectivity brain networks
% from resting state fMRI data. Medical Image Analysis, 2013,17(3):365-374 .


%% Specify input parameters
% ------------------------------------------------------------  
n_subjects = 9;

% Input information 
fs = 1; % sampling_frequency = 1/TR
temporal_mask = []; %to do motion scrubbing - already taken care of
    
% Input basis function structure
xBF.name = 'Canonical HRF (with time and dispersion derivatives)';

xBF.len = 32;                            % length in seconds of basis
xBF.order = 1;                           % order of basis set

xBF.T = 16;                              % number of subdivisions of TR
                                         % microtime resolution 
xBF.TR = 1/fs;                           % scan repetition time (seconds)
xBF.T0 = fix(xBF.T/2);                   % first time bin 
xBF.dt = xBF.TR/xBF.T;                   % length of time bin (seconds)
                                         % xBF.dt = xBF.TR/xBF.T
                                         
xBF.AR_lag = 1;                          % AR(1) noise autocorrelation
xBF.thr = 1;                             % threshold for BOLD events 

% if xBF.TR<=2; localK = 1; else; localK = 2; end
 
min_onset_search = 1;  % minimum delay allowed between event and HRF onset (seconds)
max_onset_search = 8; % maximum delay allowed between event and HRF onset (seconds)

xBF.lag = fix(min_onset_search/xBF.dt)... % array of acceptable lags (bins)
    :fix(max_onset_search/xBF.dt);

%% Compute HRF delay (seconds) and events time points for each subject 

for subj = 1:n_subjects
    
    %% Specify and load input data 
    % ------------------------------------------------------------ 
    
    % Input BOLD data 
    %path_data_in = ['/home/francisca/Documents/Tese/Implementation/Dataset/source_reconstructed_FC/fmri_connect_desikan/subj0' num2str(subj) '-7T'];
    path_data_in = ['/home/francisca/Documents/Tese/Implementation/Dataset/source_reconstructed_FC/fmri_connect_destrieux/subj0' num2str(subj) '-7T'];
    load([path_data_in '/timeseries_regr_wgm_globmean_filt0009008.mat'])
    
    %BOLD = timeseries(19:86,:); %desikan atlas
    BOLD = timeseries(13:160,:); %destrieux atlas
    
    % Load data referring the good and bad indices, identified in EEG, due to motion:
    load('/home/francisca/Documents/Tese/Implementation/Dataset/source_reconstructed_FC/eegepochs_brainstorm_manual_rejected.mat')

    ind_remove = subject(subj).sess.bad; %array with indices to remove from BOLD time series
    num_time_points = size(BOLD,2);
    num_regions = size(BOLD,1);

    aux = true(1,num_time_points); %auxiliar array to store the timepoints to keep

    for i = 1:size(ind_remove,2)
        aux(ind_remove(i)) = false;
    end
        
    data = BOLD(:,aux);
    num_time_points = length(data);
    
    %% Estimate BOLD events and HRF 
    % ------------------------------------------------------------

    % Estimate the HRF basis function and BOLD events 
    % bf - orthogonalized HRF basis function 
    % event_bold - time points of estimated pseudo-events for each region of interest
    [beta_hrf, bf, event_bold(subj).time_points] = ...
            rsHRF_estimation_temporal_basis(data',xBF,temporal_mask);

    % Scale HRF basis function - by the scaling parameter
    hrfa = bf*beta_hrf(1:size(bf,2),:);

    % Estimate HRF parameters from the 
    % hrf_pars(1) - height
    % hrf_pars(2) - time to peak (best lag) in seconds
    % hrf_pars(3) - width 
    hrf_pars = zeros(3,num_regions);
    %events_plot = zeros(num_regions, num_time_points);

    for region_id=1:num_regions
        
        hrf1 = hrfa(:,region_id);
        [hrf_pars(:,region_id)] = wgr_get_parameters(hrf1,xBF.TR/xBF.T);
        
        %events_time_points = event_bold(subj).time_points(region_id);
        %events_plot(region_id, events_time_points{1,:}) = 5;

    end
    
    %for region_id=1:num_regions
  
    %    figure
    %    plot(data(region_id,:));
    %    hold on
    %    plot(events_plot(region_id,:));
  
    %end
    
    hrf_delay(subj).all = hrf_pars(2,:); 
    hrf_delay(subj).mean = mean(hrf_pars(2,:));
    hrf_delay(subj).max = max(hrf_pars(2,:)); 
    hrf_delay(subj).min = min(hrf_pars(2,:));
    
end

save('/home/francisca/Documents/Tese/Implementation/HRF_delays_dstrx.mat', 'hrf_delay');
save('/home/francisca/Documents/Tese/Implementation/event_bold_time_points_dstrx.mat', 'event_bold');


%% Perform deconvolution of the BOLD signal
% ------------------------------------------------------------ 

%disp('Performing deconvolution ...');
%tic

%hrfa_TR = resample(hrfa,1,xBF.T);
%hrf=hrfa_TR;

% Deconvolution using a Wiener (restoration) filter
%H=fft([hrf; zeros(number_time_points-length(hrf),1)]);
%M=fft(data);
%data_deconv = ifft(conj(H).*M./(H.*conj(H)+.1*mean(H.*conj(H))));

%toc
%disp('Done');


%% Plot and save results 
% ------------------------------------------------------------  

%event_plot=sparse(1,number_time_points);
%event_plot(event_bold{1,1})=1;
%figure(1);hold on;plot((1:length(hrfa(:,1)))*xBF.TR/xBF.T,hrfa(:,1),'b');xlabel('time (s)')
%title('HRF')
%figure;
%plot((1:number_time_points)*xBF.TR/xBF.T,zscore(data)); 
%hold on;
%plot((1:n_pnts)*xBF.TR/xBF.T,zscore(data_deconv),'r');
%stem((1:n_pnts)*xBF.TR/xBF.T,event_plot,'k');
%legend('BOLD','deconvolved','events');xlabel('time (s)')