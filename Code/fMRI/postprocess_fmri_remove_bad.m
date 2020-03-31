function postprocess_fmri_remove_bad(mode)

% This function executes a post-processing step in the connectivity
% matrices so as to remove the columns and rows associated to bad indices,
% already removed from the corresponding EEG data (due to he

% with mode = 'dsk', this function uses the Desikan atlas to define the 
% regions of interest; with mode = 'dstrx', it uses the Destrieux atlas.
% This way it's already prepared for multiple options of parcellation, for
% the same post-processing

%% Removal of bad indices/time-points of the fMRI connectivity data for each subject
    
    % Load data referring the good and bad indices, identified in EEG:
    load('/home/francisca/Documents/Tese/Implementation/Datasets/source_reconstructed_FC/eegepochs_brainstorm_manual_rejected.mat')
    
    % Load fMRI connectivity data according to parcellation choice:
    if strcmp(mode,'dsk')
        % In the case parcellation chosen is according to Desikan atlas
        Directory  = '/home/francisca/Documents/Tese/Implementation/Datasets/source_reconstructed_FC/fmri_connect_desikan/';
        fMRI_data  = load('/home/francisca/Documents/Tese/Implementation/Datasets/source_reconstructed_FC/fmri_connect_desikan/fmri_results_Hilbert.mat');
    elseif strcmp(mode,'dstrx')
        % In the case parcellation chosen is according to Destrieux atlas
        fMRI_data  = load('/home/francisca/Documents/Tese/Implementation/Datasets/source_reconstructed_FC/fmri_connect_destrieux/fmri_results_Hilbert.mat');
        Directory  = '/home/francisca/Documents/Tese/Implementation/Datasets/source_reconstructed_FC/fmri_connect_destrieux/';
    end
    
    % Remove connectivity matrices for bad time points
    for s = 1:(fMRI_data.N_subjects)
    
        ind_remove = subject(s).sess.bad; %array with indices to remove from fMRI connectivity data
        aux = true(1,fMRI_data.Tmax); %auxiliar array to store the timepoints to keep

        for i = 1:size(ind_remove,2)
           aux(ind_remove(i)) = false;
        end
        
        connFMRI = fMRI_data.iFC_all(:,:,aux,s); %connectivity matrices for each good timepoint
        
        if strcmp(mode,'dsk')
            % In the case parcellation chosen is according to Desikan atlas
            save([Directory 'subj0' num2str(s) '-7T/conn_desi_phase_coh_time_fmri.mat'],'connFMRI')
        elseif strcmp(mode,'dstrx')
            % In the case parcellation chosen is according to Destrieux atlas
            save([Directory 'subj0' num2str(s) '-7T/conn_destrx_phase_coh_time_fmri.mat'],'connFMRI')
        end
       
    end

end