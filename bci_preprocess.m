classdef bci_preprocess
    %BCI_PREPROCESS Preprocessing phase of BCI
    %   This class contains all the funcions necessary in the preprocessing
    %   phase of a GDF file used in BCI
    
    properties
        % file and paths
        filename    % name of the file to be analysed
        tmpraw      % file where the relevant parameters of the GDF file are stored
        tmppsd      % file where the PSD of the buffered file is stored
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
        WinPeriod   % width of a buffering window [s]
        WinLength   % length of a window [samples]
        WinStep     % distance between following windows [samples]
        WinStart    % begin of each window [samples]
        WinStop     % end of each window [samples]
        NumWins     % number of windows in current data
        WinPos      % window of position of the events (each event is associated
                    % with the window starting when the event takes place)
        
        % Power Spectral Density parameters
        NFFT        % number of discrete Fourier transform points
        LenF        % number of frequencies
        func        % function to be used for FFT ('periodogram' or 'pwelch')
        
        % Power Spectral Density
        PSD     % power spectral density [window x frequency x channel]
        f       % frequencies of the PSD
    end
    
    methods
        function this = bci_preprocess(fname,datapath,codepath)
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

            % set buffering parameters
            this = this.set_buffering(0.5,32);
            
            % compute PSD
            this = this.compute_PSD();
        end
                
        function this = load_file(this)
            if ~exist([this.tmppath this.tmpraw],'file')            
                %load data and header
                [ls,lh] = sload([this.datapath this.filename]);
   
                % save raw data to temporary file
                save([this.tmppath this.tmpraw],'ls','lh');
            else
                tmp = load([this.tmppath this.tmpraw]);
                ls = tmp.ls;
                lh = tmp.lh;
            end
            
            this.h = lh;
            this.s = ls;

            % load event type and position
            this.type = this.h.EVENT.TYP;
            this.pos = this.h.EVENT.POS;

            % load sample rate
            this.SampleRate = this.h.SampleRate;

            % use only first 16 channels (reference is unused)
            this.s = this.s(:,1:16);
            this.NumChannels = size(this.s,2);
        end
        
        function this = set_buffering(this,period,step)
            this.WinPeriod = period;
            this.WinLength = this.WinPeriod*this.SampleRate;
            this.WinStep = step;
            this.WinStart = 1:this.WinStep:size(this.slap,1)-this.WinLength;
            this.WinStop = this.WinStart + this.WinLength;
            this.NumWins = length(this.WinStart);
            % position (window) of the events
            this.WinPos = floor(this.pos./this.WinStep);
        end
        
        function this = compute_PSD(this,func)
            if nargin < 2
                func = 'periodogram';
            end
                        
            this.NFFT = max(256,2^nextpow2(this.WinLength));
            if mod(this.NFFT,2) == 0
                this.LenF = this.NFFT/2+1;
            else
                this.LenF = (this.NFFT+1)/2;
            end
            
            if strcmp(func,'pwelch')
                this.tmppsd = [this.filename(1:end-4) '.pwelch.mat'];
            else
                this.tmppsd = [this.filename(1:end-4) '.period.mat'];
            end
            
            if ~exist([this.tmppath this.tmppsd],'file')
                if strcmp(func,'pwelch')
                    % use pwelch algorithm for PSD
                    lPSD = zeros(this.NumWins,this.LenF,this.NumChannels);

                    for win = 1:this.NumWins

                        if ~mod(win,1000)
                            disp(win)
                        end

                        % extract signal
                        sloc = this.slap(this.WinStart(win):this.WinStop(win),:);

                        % for each channel, calculate PSD
                        for ch = 1:this.NumChannels
                            [psd, lf] = pwelch(sloc(:,ch),[],[],this.NFFT,this.SampleRate);
                            lPSD(win,:,ch) = psd;
                        end
                    end
                else
                    % use periodogram for PSD
                    lPSD = zeros(this.NumWins,this.LenF,this.NumChannels);
                    
                    for win = 1:this.NumWins

                        if ~mod(win,1000)
                            disp(win)
                        end

                        % extract signal
                        sloc = this.slap(this.WinStart(win):this.WinStop(win),:);
                        
                        % generate periodogram window
                        window = hann(size(sloc,1));

                        % for each channel, calculate PSD
                        for ch = 1:this.NumChannels
                            [psd, lf] = periodogram(sloc(:,ch),window,this.NFFT,this.SampleRate);
                            lPSD(win,:,ch) = psd;
                        end
                    end
                end
                save([this.tmppath this.tmppsd],'lPSD','lf');
            else
                tmp = load([this.tmppath this.tmppsd],'lPSD','lf');
                lPSD = tmp.lPSD;
                lf = tmp.lf;
            end
            
            this.PSD = lPSD;
            this.f = lf;
            
        end


    end
    
end

