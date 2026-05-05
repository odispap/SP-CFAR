%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CFAR Detector
% For the detection of ships in SAR images
% Basic version - includes CA, Weibull, etc 
%
% September 2018
% Odysseas Pappas
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%
% CLEAN UP
%%%%
clear all;
close all;

%%%%
% INITIALISATIONS
% SET WINDOW SIZE FOR CFAR
%%%%
back_size = 30; % background band size (RADIAL)
guard_size = 10; % guard band size (RADIAL)
PFA = 0.001; % desired Probability of False Alarm
detections = 0; % initialise detection count

%%%%
% READ IMAGE
%%%%
[filename, pathname] = uigetfile({'*.jpg;*.tif;*.png;*.bmp;*.gif;*.tiff','All Image Files';...
    '*.*','All Files' },'Select images', 'MultiSelect', 'off');
% input image
img = imread([pathname,filename]);
img = double(img);
%img = im2double(imread([pathname,filename]));
figure, imshow(img, []);
img_ext = wextend(2,'sym',img,back_size); %Normal Extended Image (padded for window operations)
[h, l] = size(img);


%%%%
% APPLY THROUGH IMAGE
% This is here implemented via a loop, could be done more efficiently
% The loop contains a number of CFAR modes, comment in/out as desired
%%%%
for i = 1:h 
    for j = 1:l
        ii = i + back_size;
        jj = j + back_size;
        cut = img(i,j); % CUT is Cell Under Test (pixel of interest)
        % form bands
        guardband = img_ext(ii-guard_size:ii+guard_size, jj-guard_size:jj+guard_size);
        background = img_ext(ii-back_size:ii+back_size, jj-back_size:jj+back_size);
        backgroundNAN = background;
        backgroundNAN((back_size-guard_size):(back_size+guard_size), (back_size-guard_size):(back_size+guard_size))=NaN; %Remove CUT and GUARD
        re_back = reshape(backgroundNAN,[(2*back_size+1)*(2*back_size+1),1]); %Reshape into vector for calc
        re_back(isnan(re_back)) = [];
        
        %%%%
        % WEIBULL CFAR
        %%%%    
        pd = mle((re_back+1),'distribution','weibull'); %+1 for MLE as it doesn't like 0 intensity values
        scale = pd(1,1);
        shape = pd(1,2);
        term = nthroot(-log(PFA),shape);
        TH = scale*term; 
       
       % Apply threshold
        if cut >= TH
            maskimg(i,j) = 1;
            detections = detections + 1;
        else
            maskimg(i,j) = 0;
        end
        
    end
end


figure,imshow(maskimg,[]);
