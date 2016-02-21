classdef bci
    %BCI This class encapsulates the functions necessary for BCI
    %   This class contains all the functions for the analysis of a single
    %   GDF file for Brain Computer Interface
    
    properties (Constant)
        % events
        EVENT_TRIAL_START = 1;
        EVENT_LEFT_HAND   = 769;
        EVENT_RIGHT_HAND  = 770;
        EVENT_BOTH_FEET   = 771;
        EVENT_REST        = 783;
        EVENT_CFEEDBACK   = 781;
        TRIAL_PERIOD      = 4;    % in seconds:
                                  % - 1 seconds -> CUE 
                                  % - 3 seconds -> CONTINUOUS FEEDBACK
    end
    
    properties
        % file and paths
        filename    % name of the file to be analysed
        tmpraw     % file where the relevant parameters of the GDF file are stored
        datapath    % directory containing the data to be analysed
        codepath    % directory containing the code
        tmppath     % directory where elaboration files are saved (they can be cancelled)
        
        % EEG data from the GDF file
        s           % raw data [sample x channel]
        slap        % data after application of laplacian filter
        h           % header of the GDF file
        pos         % array containing the position (in sample) of the events
        type        % array containing the type of event
        SampleRate  % sample rate of the data [Hz]
        NumChannels % number of EEG channels
        
        % buffering parameters
        
        
    end
    
    methods
        function this = bci(fname,datapath,codepath)
          
            this.filename = fname;
            this.datapath = datapath;
            this.codepath = codepath;
            this.tmppath = [datapath 'tmp/'];
            if ~exist(this.tmppath,'dir')
                mkdir(this.tmppath);
            end
            this.tmpraw = [this.filename(1:end-4) '.raw.mat'];
            
            load([codepath 'extra/lapmask_16ch.mat'], 'lapmask');

            this = this.load_file();
               
            % apply laplacian filter
            this.slap = this.s*lapmask;

        end
                
        function this = load_file(this)
            if ~exist([this.tmppath this.tmpraw],'file')            
                %load data and header
                [s,h] = sload([this.datapath this.filename]);
   
                % save raw data to temporary file
                save([this.tmppath this.tmpraw],'s','h');
            else
                tmp = load([this.tmppath this.tmpraw]);
                s = tmp.s;
                h = tmp.h;
            end
            
            this.h = h;
            this.s = s;

            % load event type and position
            this.type = this.h.EVENT.TYP;
            this.pos = this.h.EVENT.POS;

            % load sample rate
            this.SampleRate = this.h.SampleRate;

            % use only first 16 channels (reference is unused)
            this.s = this.s(:,1:16);
            this.NumChannels = size(this.s,2);
        end
        
%         function buffering(
            

    end
    
end

