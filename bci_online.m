classdef bci_online < bci_stack
    %BCI_ONLINE Online BCI pipeline
    %   Online analysis of BCI data
    
    properties (Constant)
        EVENT_TRIAL_START = 1;
        EVENT_LEFT_HAND   = 769;
        EVENT_RIGHT_HAND  = 770;
        EVENT_BOTH_FEET   = 771;
        EVENT_REST        = 783;
        EVENT_CFEEDBACK   = 781;
        EVENT_TARGET_HIT  = 897;
        EVENT_TARGET_MISS = 898;
    end
    
    properties
        % useful data for the online cycle
        Events          % events used in the online cycle
        WinCuePos       % position of the cues
        WinEndPos       % position of the end of the trial (hit or miss)
        TrialLb         % label of the trial
        TrialResult     % result of the trial (hit or miss)        
        NumTrials       % number of trials
        
        % feature extraction
        f_PSD           % frequencies used as features
        Lenf_PSD        % number of frequencies used as features
        logPSD_sel      % PSD for the selected features
        NumObs          % number of observations
        NumFeatures     % number of features        
        
        % posterior probabilities
        selFeatures     % feature selected for the classification (from offline analysis)
        class           % classifier used for the BCI
        pp              % posterior probability for the window [window x pp]

        % accumulation framework parameters
        alpha           % alpha used in the accumulation framework
        P               % accumulated probability for the trial [trial]
        threshold       % threshold used in the accumulation framework
        result          % result of the accumulation framework
        
        % results of the online cycle
        hit             % events where cue and result are equal
        miss            % events where cue and result are different
        none            % events where no result is issued
        probSucc        % success probability of the cycle
    end
    
    methods
        
        function this = bci_online(fnames,datapath,codepath)
            this@bci_stack(fnames,datapath,codepath);
            
            this.Events = [this.EVENT_LEFT_HAND this.EVENT_BOTH_FEET];
            
            this.WinCuePos = this.WinPos(this.type==this.Events(1) | this.type==this.Events(2));
            this.WinEndPos = this.WinPos(this.type==this.EVENT_TARGET_HIT | this.type==this.EVENT_TARGET_MISS);
            
            this.NumTrials = length(this.WinCuePos);
            this.TrialLb = this.type(this.type==this.Events(1) | this.type==this.Events(2));
            this.TrialResult = this.type(this.type==this.EVENT_TARGET_HIT | this.type==this.EVENT_TARGET_MISS);
        end
        
        function this = set_classifier(this,cls)
            this.class = cls;
        end
        
        function this = featureExtraction(this,freq,selFeat)
            [this.f_PSD,indPSD,~] = intersect(this.f,freq);
            this.Lenf_PSD = length(this.f_PSD);
            
            PSD_tmp = this.PSD(:,indPSD,:);
            
            this.NumFeatures = this.Lenf_PSD * this.NumChannels;
            this.NumObs = this.NumWins;
            
            PSD_tmp = permute(PSD_tmp, [1 4 2 3]);
            PSD_tmp = reshape(PSD_tmp, [this.NumObs,this.NumFeatures]);
            PSD_tmp = log(1 + PSD_tmp);
            
            this.selFeatures = selFeat;
            this.logPSD_sel = PSD_tmp(:,this.selFeatures);
        end
        
        function this = computePP(this,cls)
            addpath(genpath([this.codepath '/extra/classification']));
            
            this = this.set_classifier(cls);
            [this.pp,cls] = classTest(this.class,this.logPSD_sel);
        end
        
        function this = accumParam(this,a,thres)
            this.alpha = a;
            this.threshold = thres;
        end
        
        function this = onlineCycle(this)
            this.P = zeros(this.NumTrials,2);
            this.result = zeros(this.NumTrials,1);
            for i = 1:length(this.WinCuePos)
%                 disp(['Cycle number ' num2str(i) ', Event = ',num2str(this.TrialLb(i))]);
                this.P(i,:) = this.pp(this.WinCuePos(i),:);
%                 disp(['P1 = ' num2str(this.P(i,1)) ' - P2 = ' num2str(this.P(i,2))]);
                for j = this.WinCuePos(i)+1:this.WinEndPos(i)
                    if this.P(i,1)>this.threshold || this.P(i,2)>this.threshold
                        break
                    end
                    this.P(i,1) = this.alpha*this.P(i,1) + (1-this.alpha)*this.pp(j,1);
                    this.P(i,2) = this.alpha*this.P(i,2) + (1-this.alpha)*this.pp(j,2);
%                     disp(['P1 = ' num2str(this.P(i,1)) ' - P2 = ' num2str(this.P(i,2))]);
                end
                if this.P(i,1)>this.threshold
%                     disp(['Event = ' num2str(this.Events(1))]);
                    this.result(i) = this.Events(1);
                elseif this.P(i,2)>this.threshold
%                     disp(['Event = ' num2str(this.Events(2))]);
                    this.result(i) = this.Events(2);
                else
%                     disp(['Impossible to discriminate'])
                end
            end
            
            this.hit = length(find(this.result==this.TrialLb));
            this.none = length(find(this.result==0));
            this.miss = this.NumTrials - this.none - this.hit;
            
            this.probSucc = this.hit/this.NumTrials;
            
        end

    end
    
end

