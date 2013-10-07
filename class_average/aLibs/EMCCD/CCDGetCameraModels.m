function s=CCDGetCameraModels
% function models=CCDGetCameraModels
% Pick up the models structure describing our cameras.



% First, find our local directory
pa=fileparts(which('CCDGetCameraModels'));
% Retrieve parameters from SpectModel.mat in the local directory
if numel(pa)<1
    pa='.';
end;
load([pa '/SpectModel.mat']);
s=models;