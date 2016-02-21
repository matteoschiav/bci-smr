
rootpath = '/home/schiavon/Matteo/luxor/cisas/bci/';
trainpath = [rootpath 'bci4neurorobotics-data/smr/20150914_b1/'];
codepath = [rootpath 'bci4neurorobotics-code/'];

% oldfolder = cd([codepath 'biosig']);
% biosig_installer;
% cd(oldfolder);

filetrain{1} = 'b1.20150914.102659.offline.mi.mi_rlbf.gdf';

offline = bci(filetrain{1},trainpath,codepath);
offline = offline.compute_PSD();