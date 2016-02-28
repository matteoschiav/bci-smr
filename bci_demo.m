clc
clear all
close all
clear classes

rootpath = '/home/schiavon/Matteo/luxor/cisas/bci/';
trainpath = [rootpath 'bci4neurorobotics-data/smr/20150914_b1/'];
codepath = [rootpath 'bci4neurorobotics-code/'];

% oldfolder = cd([codepath 'biosig']);
% biosig_installer;
% cd(oldfolder);

% classifier train
filetrain{1} = 'b1.20150914.102659.offline.mi.mi_rlbf.gdf';
filetrain{2} = 'b1.20150914.104000.offline.mi.mi_rlbf.gdf';

% classifier test
filetrain{3} = 'b1.20150914.105054.offline.mi.mi_rlbf.gdf';

% TEST DATA
% trainpath = [rootpath 'bci4neurorobotics-data/smr/20150831_b4/'];
% 
% filetrain{1} = 'b4.20150831.121033.offline.mi.mi_rlsf.gdf';
% filetrain{2} = 'b4.20150831.123916.offline.mi.mi_rlsf.gdf';
% filetrain{3} = 'b4.20150831.125323.offline.mi.mi_rlsf.gdf';
% filetrain{4} = 'b4.20150831.132126.offline.mi.mi_rlsf.gdf';
% END TEST DATA


offline = bci_offline(filetrain,trainpath,codepath);
offline = offline.calcERD(4:48);
% offline.plotERD();
% offline.plotPSD();

%%
offline = offline.paramFE(4:48,[2 5],{'left hand' 'both feet'});
offline = offline.divideTrainTest(60);
offline = offline.computeDP();
offline.plotDP();
offline.plotCanonical();
offline = offline.feature_extraction(4);

classifiers = {'lda','qda','knn','gau'};
class_desc = {'Linear Discriminant Analysis',...
    'Quadratic Discriminant Analysis',...
    'k-nearest neighbour Classification',...
    'Gaussian Classification'};
offclass = {};

for i = 1:length(classifiers)
    tmpclass = offline.trainClass(classifiers{i});
    offclass{i} = tmpclass;
    disp(class_desc{i});
    fprintf('Success probability on train data: %f\n',offclass{i}.probSuccTrain);
    offclass{i} = offclass{i}.testClass();
    fprintf('Success probability on test data: %f\n',offclass{i}.probSuccTest);
end


%% online analysis
onlinepath = [rootpath 'bci4neurorobotics-data/smr/20150915_b1/'];

onfile{1} = 'b1.20150915.101044.online.mi.mi_rlbf.gdf';
onfile{2} = 'b1.20150915.101734.online.mi.mi_rlbf.gdf';
onfile{3} = 'b1.20150915.102331.online.mi.mi_rlbf.gdf';
onfile{4} = 'b1.20150915.103128.online.mi.mi_rlbf.gdf';

online = {};
for i = 1:length(classifiers)
    online{i} = bci_online(onfile,onlinepath,codepath);
    online{i} = online{i}.featureExtraction(offclass{i}.f_FE,offclass{i}.selFeatures);
    online{i} = online{i}.computePP(offclass{i}.class);
end
