%% Initialization

% Folder where data is saved.
Directory  = 'D:\Tese\EEG\';

% Code to run EEGLAB toolbox
% global main_path; main_path = '/home/aaires/Toolbox/';
% addpath(genpath([ main_path, 'EEGLAB' ]))

% Data parameters
fmri_TR = 1;
N_subjects = 9;
Tmax = 470;
nb_channels = 59;

psd_all = zeros(2^12/2+1,N_subjects*Tmax,nb_channels);
bad_TR  = cell(1,N_subjects);

for s = 1:N_subjects
    
    %Load EEG data
    load([Directory 'EEG' num2str(s) '.mat'])
    
    
    TR = find(eeg.Mt(:,1));     % eeg timepoint that marks each TR
    TR = TR(eeg.dummies+1:end); % removal of 10 first TR (equivalent to what was done in fMRI)
    TR = [TR;size(eeg.Yb,1)];
    
    
    % Identification of bad segments
    
    bad_points = find(eeg.Mt(:,5)==1);
    for i = 1:length(bad_points)
        ok = false;
        j = 1;
        while ~ok
            if (bad_points(i) - TR(j)) < 0
                ok = true;
            else
                j = j+1;
            end
        end
        bad_TR{s}.idx(i) = j;
    end
    
    % computation of psd for each TR using 2^12 points in fft
       
    psd_eeg=zeros(2^12/2+1,Tmax,nb_channels);
    
    for i=1:nb_channels
        
        for j=1:Tmax
            fft_eeg=fft(eeg.Yb(TR(j):TR(j+1)-1,i),2^12);
            N=length(eeg.Yb(TR(j):TR(j+1)-1));
            fft_eeg=fft_eeg(1:2^12/2+1);
            psd_eeg_ch=(1/500)*abs(fft_eeg).^2;
            psd_eeg_ch(2:end-1) = 2 * psd_eeg_ch(2:end-1);
            psd_eeg(:,j,i)=psd_eeg_ch;
        end        
    end
    
    psd_all(:,(s-1)*Tmax+1:Tmax*s,:)=psd_eeg;
    
    disp(['Subj ' num2str(s) ' completed'])
end

% Useful parameters:
fs = eeg.Fs;
chanlocs = eeg.chanlocs;
freq = 0:fs/2^12:fs/2;

%% Compute problematic TRs

bad_TRs = [];
for s = 1:N_subjects
%     aux = unique(bad_TR{s}.idx);
%     l = zeros(1,length(aux));
%     for i = 1:length(aux)
%         l(i) = length(find(bad_TR{s}.idx == aux(i)));
%     end
%     aux(l<500/2) = [];
    bad_TRs = [bad_TRs 470*(s-1)+unique(bad_TR{s}.idx)];
end

% Compute good_TRs
t_all = 1:Tmax*N_subjects;
good_TRs = t_all(~ismember(t_all,bad_TRs));

save('eeg_analysis.mat','psd_all','bad_TRs','chanlocs','fmri_TR','freq','fs','good_TRs','N_subjects','nb_channels','Tmax','t_all','-v7.3')