function fmri_phase_coherence(mode)

% This function processes the BOLD data so as to obtain the phase coherence
% as a connectivity metric
% (code adapted from Joana Cabral)

% 1 - Read the BOLD data from the folders and computes the BOLD phases
% 2 - Calculate the instantaneous BOLD synchronization matrix

% with mode = 'dsk', this function uses the Desikan atlas to define the 
% regions of interest; with mode = 'dstrx', it uses the Destrieux atlas.

 %% Initialization
    
    % Data parameters:
    N_subjects = 9;
    fmri_TR    = 1;
    Tmax       = 480;
   
    % Folder where data is saved
    if strcmp(mode,'dsk')
        % In the case parcellation chosen is according to Desikan atlas
        N_areas    = 68;
        Directory  = '/home/francisca/Documents/Tese/Implementation/Datasets/source_reconstructed_FC/fmri_connect_desikan/';
    elseif strcmp(mode,'dstrx')
        % In the case parcellation chosen is according to Destrieux atlas
        N_areas    = 148;
        Directory  = '/home/francisca/Documents/Tese/Implementation/Datasets/source_reconstructed_FC/fmri_connect_destrieux/';
    end
            
    iFC_all  = zeros(N_areas,N_areas,Tmax,N_subjects); % 4D matrix to store PC for each TR and subject
    
    for s = 1:N_subjects
        % USER: Adapt to load here the BOLD time-series matrix (NxT) from each subject
        load([Directory 'subj0' num2str(s) '-7T/timeseries_regr_wgm_globmean_filt0009008.mat'])
        % using regressed and bandpass filtered time series

        % Remove extra ROIs from BOLD time series: 
        if strcmp(mode,'dsk')      
            % (first 18 - Desikan atlas)
            BOLD = timeseries(19:86,:); % N_areas x Tmax
        elseif strcmp(mode,'dstrx')
            % (first 12 - Destrieux atlas)
            BOLD = timeseries(13:160,:); % N_areas x Tmax
        end 
        
        % Get the BOLD phase using the Hilbert transform:
       
        Phase_BOLD = zeros(N_areas,Tmax);
        
        for seed = 1:N_areas
            BOLD(seed,:) = BOLD(seed,:)-mean(BOLD(seed,:));
            Phase_BOLD(seed,:) = angle(hilbert(BOLD(seed,:)));
        end
        
        % Computing Phase Coherence for all TRs and subjects
        
        for t = 1:Tmax
            
            % Calculate the Instantaneous FC (BOLD Phase Synchrony) between
            % any two regions
            iFC = zeros(N_areas, N_areas);
            for n = 1:N_areas
                for p = 1:N_areas
                    iFC(n,p) = cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
                end
            end
            iFC_all(:,:,t,s) = iFC;
        end  
    end
    save([Directory 'fmri_results_Hilbert.mat'],'N_subjects', 'N_areas', 'fmri_TR', 'Tmax','iFC_all')
end