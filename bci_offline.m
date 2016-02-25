classdef bci_offline
    %BCI_OFFLINE Offline BCI pipeline
    %   Offline analysis of the BCI data: trial extraction, feature
    %   selection, feature extraction and classifier training.
    
    properties (Constant)
        EVENT_TRIAL_START = 1;
        EVENT_LEFT_HAND   = 769;
        EVENT_RIGHT_HAND  = 770;
        EVENT_BOTH_FEET   = 771;
        EVENT_REST        = 783;
        EVENT_CFEEDBACK   = 781;
    end
    
    properties
        % data for the offlne analysis
        data        % class bci_stack containing the data
        
        % files and directories
        datapath    % directory where data are stored
        codepath    % directory where code is stored
        filenames   % name of files belonging to the dataset
        
        % trial parameters
        TRIAL_START     % start of the trial (from cue) [s]
        TRIAL_STOP      % stop of the trial (from cue) [s]
        TRIAL_PERIOD    % trial duration [s]
        
        % trial selection
        TrialId         % identifier of the selected trials
        TrialLb         % label of each event
        CuePos          % position of each event [sample]
        CueWin          % position of each event [window]
        NumTrials       % number of the selected trials
        TrialStartWin   % begin of the trial [window]
        TrialStopWin    % end of the trial [window]
        WinPerTrial     % trial duration [window]
        CueInTrial      % position of the cue in the trial [window]
        
        % trial PSD
        PSD_trial       % PSD divided by trials [window x frequency x channel x trial]
        PSD_baseline    % PSD of the baseline (from trial start to cue)
        PSD_task        % PSD of the task (from cue to end of trial)
        
        % ERD|ERS data
        ERD             % ERD|ERS
        f_ERD            % frequencies where ERD|ERS is computed
    end
    
    methods
        
        function this = bci_offline(fnames,datapath,codepath)
            this.filenames = fnames;
            this.datapath = datapath;
            this.codepath = codepath;
            
            this.data = bci_stack(fnames,datapath,codepath);
            
            % set trial parameters
            this = this.trial_parameters(-1,4,1);
            
            % extract trials
            sel_trials = {'left hand' 'right hand' 'both feet'};
            this = this.trial_extraction(sel_trials);
            
        end
        
        function this = trial_parameters(this,start,stop,cueintrial)
            this.TRIAL_START = start;
            this.TRIAL_STOP = stop;
            this.TRIAL_PERIOD = this.TRIAL_STOP - this.TRIAL_START;
            this.WinPerTrial = this.TRIAL_PERIOD.*this.data.WinRate;
            this.CueInTrial = cueintrial.*this.data.WinRate;
        end
        
        function this = trial_extraction(this, sel_trials)
            % Trials that can be selected:
            %  * left hand
            %  * right hand
            %  * both feet
            %  * rest
            if nargin < 2
                sel_trials = {'left hand' 'right hand' 'both feet' 'rest'};
            end
            
            id_num = [this.EVENT_LEFT_HAND this.EVENT_RIGHT_HAND ...
                this.EVENT_BOTH_FEET this.EVENT_REST];
            id_string = {'left hand' 'right hand' 'both feet' 'rest'};
            
            % extract selected trials
            this.TrialId = zeros(size(this.data.type));
            for i = 1:length(id_num)
                if ~isempty(ismember(sel_trials,id_string(i)))
                    this.TrialId = this.TrialId | (this.data.type == id_num(i));
                end
            end
            
            this.TrialLb = this.data.type(this.TrialId);
            this.CuePos = this.data.pos(this.TrialId);
            this.CueWin = floor(this.CuePos./this.data.WinStep);
            this.NumTrials = length(this.CueWin);
            
            this.TrialStartWin = this.CueWin + this.TRIAL_START.*this.data.WinRate;
            this.TrialStopWin = this.CueWin + this.TRIAL_STOP.*this.data.WinRate - 1;
            
            % extract trials
            this.PSD_trial = zeros(this.WinPerTrial,this.data.LenF,this.data.NumChannels,this.NumTrials);
            for i = 1:this.NumTrials
                this.PSD_trial(:,:,:,i) = ...
                    this.data.PSD(this.TrialStartWin(i):this.TrialStopWin(i),:,:);
            end
            
            % extract baseline and task from trial
            this.PSD_baseline = this.PSD_trial(1:this.CueInTrial,:,:,:);
            this.PSD_task = this.PSD_trial(this.CueInTrial+1:end,:,:,:);
        end
        
    end
    
end

