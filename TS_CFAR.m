%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Truncated Statistics CFAR Detector
% For the detection of ships in SAR images
% 
% Based on Tao et al, TGRS 2016
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
img = im2double(imread([pathname,filename]));
figure, imshow(img, []);
img_ext = wextend(2,'sym',img,back_size); %Normal Extended Image (padded for window operations)
[h, l] = size(img);
TH = zeros(h,l);


%%%%
% APPLY THROUGH IMAGE
% This is here implemented via a loop, could be done more efficiently
% The loop contains a number of CFAR modes, comment in/out as desired
%%%%
for i = 1:h %swapped l and h around - NOTE!
    for j = 1:l
        ii = i + back_size;
        jj = j + back_size;
        cut = img(i,j);
        
        % Get Bands
        guardband = img_ext(ii-guard_size:ii+guard_size, jj-guard_size:jj+guard_size);
        background = img_ext(ii-back_size:ii+back_size, jj-back_size:jj+back_size);
        backgroundNAN = background;
        backgroundNAN((back_size-guard_size):(back_size+guard_size), (back_size-guard_size):(back_size+guard_size))=NaN; %Remove CUT and GUARD
        re_back = reshape(backgroundNAN,[(2*back_size+1)*(2*back_size+1),1]); %Reshape for calc
%         re_back = re_back + 1; %remove 0 
          
        sample_mean = nanmean(re_back(:));
        trunc_thresh = 0.75*max(max(re_back));
        trunc_sample = re_back; %Copy
        trunc_sample(trunc_sample > trunc_thresh) = 0; %Truncate
        trunc_sample_mean = nanmean(trunc_sample(:));
        m = trunc_sample_mean;
        TH = -m*log(PFA);


        
       
        % Apply Threshold
        if cut >= TH
            maskimg(i,j) = 1;
            detections = detections + 1;
        else
            maskimg(i,j) = 0;
        end;

  %toc
    end;
end;

figure,imshow(maskimg,[]);

