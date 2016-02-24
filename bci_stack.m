classdef bci_stack
    %BCI_STACK Data belonging to a single set
    %   Class containing the arrays of data belonging to a single dataset
    %   (offline/online).
    
    properties
        % classes containing the dataset
        stacks      % classes belonging to the dataset
        
        % files and directories
        datapath    % directory where data are stored
        codepath    % directory where code is stored
        filenames   % name of files belonging to the dataset
        
        % EEG data
        slap        % data after application of laplacian filter [sample x channel]
        pos         % position of the events
        type        % type of the events
        
        % EEG parameters
        SampleRate  % sample rate of the data [Hz]
        NumChannels % number of EEG channels
        
        % buffering parameters
        WinPeriod   % width of a buffering window [s]
        WinLength   % length of a window [samples]
        WinStep     % distance between following windows [samples]
        WinStart    % begin of each window [samples]
        WinStop     % end of each window [samples]
        NumWins     % number of windows in current dataset
        WinRate     % rate of the windows [Hz]
        WinPos      % window of position of the events (each event is associated
                    % with the window starting when the event takes place)
        
        % PSD data
        PSD         % power spectral density of the whole dataset [window x frequency x channel]
        f           % frequencies of the PSD
        LenF        % number of frequencies of the PSD

    end
    
    methods
        function this = bci_stack(fnames,datapath,codepath)
            this.datapath = datapath;
            this.codepath = codepath;
            this.filenames = fnames;
            
            for i = 1:length(fnames)
                this.stacks{i} = bci_preprocess(fnames{i},datapath,codepath);
            end
            
            %
            
            % load EEG parameters
            this.SampleRate = this.stacks{1}.SampleRate;
            this.NumChannels = this.stacks{1}.NumChannels;
            
            % load EEG data and buffering parameters
            for i = 1:length(this.stacks)
                this.type = cat(1,this.type,this.stacks{i}.type);
                this.pos = cat(1,this.pos,this.stacks{i}.pos+size(this.slap,1));
                
                % buffering parameters
                this.WinStart = cat(1,this.WinStart,this.stacks{i}.WinStart'+size(this.slap,1));
                this.WinStop = cat(1,this.WinStop,this.stacks{i}.WinStop'+size(this.slap,1));
                this.WinPos = cat(1,this.WinPos,this.stacks{i}.WinPos+size(this.PSD,1));
                
                % load EEG data (filtered)
                this.slap = cat(1,this.slap,this.stacks{i}.slap);
                
                % load PSD data
                this.PSD = cat(1,this.PSD,this.stacks{i}.PSD);
            end
            % other buffering parameters
            this.WinPeriod = this.stacks{1}.WinPeriod;
            this.WinLength = this.stacks{1}.WinLength;
            this.WinStep = this.stacks{1}.WinStep;
            
            this.NumWins = length(this.WinStart);
            this.WinRate = this.SampleRate./this.WinStep;
                        
            % load PSD parameters
            this.f = this.stacks{1}.f;
            this.LenF = this.stacks{1}.LenF;
        end
    end
    
end

