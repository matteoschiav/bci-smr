classdef bci_offline < bci_stack
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
        
        % trial parameters
        TRIAL_START     % start of the trial (from cue) [s]
        TRIAL_STOP      % stop of the trial (from cue) [s]
        TRIAL_PERIOD    % trial duration [s]
        
        % selected event
        SelEvents       % array containing the code of the selected events
        SelEventsLb     % array containing the label of the selected events
        
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
        ERD             % ERD|ERS [%]
        actERD          % event-processed ERD|ERS [%]
        f_ERD           % frequencies where ERD|ERS is computed
        
        % feature extraction
        classFE         % classes used for feature extraction
        classFElb       % label of classes used for feature extraction
        f_FE            % frequencies used in feature extraction
        Lenf_FE         % number of frequencies used in FE
        WinFE           % task windows for feature extraction
        NumWinFE        % number of windows used in FE
        logPSD_FE       % logarithm of the PSD used for feature extraction
        vecTrial_FE     % vector containing the selected trials
        SelTrialsFE     % trials selected for FE
        NumSelTrialsFE  % number of trials selected for FE
        NumFeatures     % number of features
        NumObs          % number of observation
        
        % discriminant power
        matDP           % discriminant power
        freqDP          % frequencies used for calculating discriminant power
        dp              % index discriminability power [%]
        v               % eigenvectors
        disc            % transformed data
        
        % classification data
        selFeatures     % features selected for classification
        logPSD_sel      % PSD for the selected features
        class           % classifier
        
        % classification test (on train data)
        probSuccTrain   % success probability on train data
    end
    
    methods
        
        function this = bci_offline(fnames,datapath,codepath)            
            this@bci_stack(fnames,datapath,codepath);
                     
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
            this.WinPerTrial = this.TRIAL_PERIOD.*this.WinRate;
            this.CueInTrial = cueintrial.*this.WinRate;
        end
        
        function this = trial_extraction(this, sel_trials)
            % Trials that can be selected:
            %  * left hand
            %  * right hand
            %  * both feet
            %  * rest
            if nargin < 1
                sel_trials = {'left hand' 'right hand' 'both feet' 'rest'};
            end
            
            id_num = [this.EVENT_LEFT_HAND this.EVENT_RIGHT_HAND ...
                this.EVENT_BOTH_FEET this.EVENT_REST];
            id_string = {'left hand' 'right hand' 'both feet' 'rest'};
            
            % extract selected trials
            this.SelEvents = [];
            this.SelEventsLb = [];
            this.TrialId = zeros(size(this.type));
            for i = 1:length(id_num)
                if ~isempty(find(ismember(sel_trials,id_string{i}),1))
                    this.TrialId = this.TrialId | (this.type == id_num(i));
                    this.SelEvents = [this.SelEvents id_num(i)];
                    this.SelEventsLb{end+1} = id_string{i};
                end
            end
            
            this.TrialLb = this.type(this.TrialId);
            this.CuePos = this.pos(this.TrialId);
            this.CueWin = floor(this.CuePos./this.WinStep);
            this.NumTrials = length(this.CueWin);
            
            this.TrialStartWin = this.CueWin + this.TRIAL_START.*this.WinRate;
            this.TrialStopWin = this.CueWin + this.TRIAL_STOP.*this.WinRate - 1;
            
            % extract trials
            this.PSD_trial = zeros(this.WinPerTrial,this.LenF,this.NumChannels,this.NumTrials);
            for i = 1:this.NumTrials
                this.PSD_trial(:,:,:,i) = ...
                    this.PSD(this.TrialStartWin(i):this.TrialStopWin(i),:,:);
            end
            
            % extract baseline and task from trial
            this.PSD_baseline = this.PSD_trial(1:this.CueInTrial,:,:,:);
            this.PSD_task = this.PSD_trial(this.CueInTrial+1:end,:,:,:);
        end
        
        function this = calcERD(this,freqs)
            [this.f_ERD,ind_f,a] = intersect(this.f,freqs);
            
            baseline = permute(repmat(squeeze(mean(log(1+this.PSD_baseline))),1,1,1,80),[4 1 2 3]);
            logPSD_trial = log(1+this.PSD_trial);
            
            this.ERD = (logPSD_trial - baseline)./baseline.*100;
            
            this.ERD = this.ERD(:,ind_f,:,:);
            
            this.actERD = zeros(size(this.ERD,1),size(this.ERD,2),this.NumChannels,length(this.SelEvents));
            for ev=1:length(this.SelEvents)
                this.actERD(:,:,:,ev) = mean(this.ERD(:,:,:,this.TrialLb==this.SelEvents(ev)),4);
            end
            
        end
        
        function this = plotERD(this)
            time = this.TRIAL_START:1/this.WinRate:this.TRIAL_STOP;
            this.plot2dElectrodes(time,this.f_ERD,this.actERD);
        end
        
        function plot2dElectrodes(this,xdata,ydata,zdata)
            NumRows = 4;
            NumCol = 5;
            
            ChanLayout = [0 0 1 0 0; 2 3 4 5 6; 7 8 9 10 11; 12 13 14 15 16];
            
            Electrode = {'Fz' 'FC3' 'FC1' 'FCz' 'FC2' 'FC4' 'C3' 'C1' 'Cz' 'C4' 'C2' 'CP3' 'CP1' 'CPz' 'CP2' 'CP4'};

            for CurEv=1:length(this.SelEvents)
                figure
                for ch=1:this.NumChannels
                    [x,y] = find(ChanLayout == ch);
                    pos = (x-1)*NumCol+y;
                    subplot(NumRows,NumCol,pos);

                    imagesc(xdata,ydata,zdata(:,:,ch,CurEv)');
                    colormap(jet);

                    xlabel('[s]');
                    ylabel('[Hz]');
                    title(Electrode(ch));

                end

                subplot(NumRows,NumCol,5);
                colorbar;
                axis off;
            end
        end
        
        function plotPSD(this)
            NumRows = 4;
            NumCol = 5;
            
            ChanLayout = [0 0 1 0 0; 2 3 4 5 6; 7 8 9 10 11; 12 13 14 15 16];
            
            Electrode = {'Fz' 'FC3' 'FC1' 'FCz' 'FC2' 'FC4' 'C3' 'C1' 'Cz' 'C4' 'C2' 'CP3' 'CP1' 'CPz' 'CP2' 'CP4'};
            
            [freqs,indf,~] = intersect(this.f,4:24);

            for CurEv=1:length(this.SelEvents)
                figure
                for ch=1:this.NumChannels
                    [x,y] = find(ChanLayout == ch);
                    pos = (x-1)*NumCol+y;
                    subplot(NumRows,NumCol,pos);

                    plot(freqs,mean(mean(20.*log10(this.PSD_baseline(:,indf,ch,this.TrialLb==this.SelEvents(CurEv))),4),1),'Linewidth',2);
                    hold on
                    plot(freqs,mean(mean(20.*log10(this.PSD_task(:,indf,ch,this.TrialLb==this.SelEvents(CurEv))),4),1),'Linewidth',1);
                    xlim([min(freqs) max(freqs)])

                    xlabel('[s]');
                    ylabel('PSD');
                    title(Electrode(ch));
                end
                subplot(NumRows,NumCol,5);
                axis off;
            end
        end
        
        function this = paramFE(this,freqs,interval,classes)
            [this.f_FE,indFE,~] = intersect(this.f,freqs);
            this.Lenf_FE = length(this.f_FE);
            
            this.WinFE = interval(1)*this.WinRate+1:interval(end)*this.WinRate;
            this.NumWinFE = length(this.WinFE);
            
            [~,indCl,~] = intersect(this.SelEventsLb,classes);
            this.classFE = this.SelEvents(indCl);
            this.classFElb = this.SelEventsLb(indCl);
            
            this.SelTrialsFE = (this.TrialLb==this.classFE(1) | this.TrialLb==this.classFE(2));
            
            % put the PSD in the form required for FE
            PSD_tmp = this.PSD_trial(this.WinFE,indFE,:,this.SelTrialsFE);
            this.NumSelTrialsFE = size(PSD_tmp,4);
            
            this.NumFeatures = this.Lenf_FE*this.NumChannels;
            this.NumObs = this.NumWinFE*this.NumSelTrialsFE;
            
            PSD_tmp = permute(PSD_tmp, [1 4 2 3]);
            PSD_tmp = reshape(PSD_tmp, [this.NumObs,this.NumFeatures]);
            this.logPSD_FE = log(1 + PSD_tmp);
            
            this.vecTrial_FE = this.TrialLb(this.SelTrialsFE);
            this.vecTrial_FE = repmat(this.vecTrial_FE',[this.NumWinFE 1]);
            this.vecTrial_FE = reshape(this.vecTrial_FE,[this.NumObs 1]);
        
        end
        
        function this = computeDP(this)
            addpath([this.codepath '/extra/cva']);
            [this.dp,~,this.v,~,this.disc] = cva_tun_opt(this.logPSD_FE,this.vecTrial_FE);
            
            this.freqDP = this.f_FE;
            this.matDP = reshape(this.dp,[this.Lenf_FE this.NumChannels]);
        end
        
        function plotDP(this)
            % feature selection matrix
            figure();
            imagesc(this.freqDP,1:this.NumChannels,this.matDP');
            colormap(jet);
            colorbar;
            
            xlabel('Frequencies [Hz]');
            ylabel('Channels');
            title([ this.classFElb{1} '-' this.classFElb{2} ]);
        end
        
        function plotCanonical(this)
            addpath([this.codepath '/extra/gkde']);
            
            x1 = this.disc(this.vecTrial_FE == this.classFE(1));
            x2 = this.disc(this.vecTrial_FE == this.classFE(2));
            
            p1 = gkdeb(x1);
            p2 = gkdeb(x2);
            
            figure();
            plot(p1.x,p1.f,p2.x,p2.f,'r');
            legend(this.classFElb{1},this.classFElb{2});
        end
        
        function this = feature_extraction(this,n_feat)
            [~,sortIndex] = sort(this.dp(:),'descend');
            this.selFeatures = sortIndex(1:n_feat);
            
            this.logPSD_sel = this.logPSD_FE(:,this.selFeatures);
        end
        
        function this = trainClass(this,type)
            addpath(genpath([this.codepath '/extra/classification']));
            
            % train classifier
            this.class = classTrain(this.logPSD_sel,this.vecTrial_FE,type);
            
            % test classifier on train data
            [~, cls] = classTest(this.class,this.logPSD_sel);
            this.probSuccTrain = length(cls(cls==this.vecTrial_FE))/size(cls,1);
        end
        
    end
    
end

