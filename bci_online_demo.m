clc
clear all
close all
clear classes

rootpath = '/home/schiavon/Matteo/luxor/cisas/bci/';
% datapath = [rootpath 'bci4neurorobotics-data/smr/20150915_b1/'];
codepath = [rootpath 'bci4neurorobotics-code/'];

% oldfolder = cd([codepath 'biosig']);
% biosig_installer;
% cd(oldfolder);

onlinepath = [rootpath 'bci4neurorobotics-data/smr/20150915_b1/'];

onfile{1} = 'b1.20150915.101044.online.mi.mi_rlbf.gdf';
onfile{2} = 'b1.20150915.101734.online.mi.mi_rlbf.gdf';
onfile{3} = 'b1.20150915.102331.online.mi.mi_rlbf.gdf';
onfile{4} = 'b1.20150915.103128.online.mi.mi_rlbf.gdf';

classpath = [onlinepath '/tmp'];

classifiers = {'lda','qda','knn','gau'};
class_desc = {'Linear Discriminant Analysis',...
    'Quadratic Discriminant Analysis',...
    'k-nearest neighbour Classification',...
    'Gaussian Classification'};

online = {};
for i = 1:length(classifiers)
    classname = ['class_' classifiers{i} '.mat'];
    
    disp(class_desc{i});
    load([classpath classname],'f_FE','selFeatures','class');
    
    online{i} = bci_online(onfile,onlinepath,codepath);
    online{i} = online{i}.featureExtraction(f_FE,selFeatures);
    online{i} = online{i}.computePP(class);
    
    online{i} = online{i}.accumParam(0.05,0.995);
    online{i} = online{i}.onlineCycle();
    
    fprintf('Success probability = %f\n', online{i}.probSucc);
end
