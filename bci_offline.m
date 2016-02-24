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
        
        % repartition of data into trials
        
    end
    
    methods
    end
    
end

