%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SP-CFAR Detector (Weibull)
% For the detection of ships in SAR images
% 
% Based on Pappas et al, GRSL 2018
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
%%%%
back_size = 30; % background band size (RADIAL)
guard_size = 10; % guard band size (RADIAL)
PFA = 0.01; % desired Probability of False Alarm
detections = 0; % initialise detection count
SLIC_k = 300; % SLIC k parameter - number of superpixels
SLIC_m = 10; % SLIC m parameter - compactness

%%%%
% READ IMAGE
%%%%
[filename, pathname] = uigetfile({'*.jpg;*.tif;*.png;*.bmp;*.gif;*.tiff','All Image Files';...
    '*.*','All Files' },'Select images', 'MultiSelect', 'off');
% Input image
img = imread([pathname,filename]);
img = double(img);
% img = im2double(imread([pathname,filename]));
figure, imshow(img, []);
img_ext = wextend(2,'sym',img,back_size); %Normal Extended Image (padded for window operations)
[h, l] = size(img);
TH = zeros(h,l);
% Create 3-channel version for standard SLIC algorithm
fakeRGBimg(:,:,1) = img;
fakeRGBimg(:,:,2) = img;
fakeRGBimg(:,:,3) = img;


%%%%
% SLIC Segmentation
%%%%
[label_img, Am, Attribute] = slic(fakeRGBimg, SLIC_k, SLIC_m, 3, 'median');
figure,imshow(drawregionboundaries(label_img, img, [255]),[0 255]);

%%%%
% SPLIT INTO SUPERPIXEL NEIGHBOURHOODS
%%%%
[~, count] = size(Attribute); %actual number of generated spixels
super{count} = [];
neighbours{count}= [];
for i = 1:l
    for j = 1:h
        % What superpixel should this pixel be part of? 
        current_count = label_img(j,i);
        % Add to the population of that superpixel
        super{current_count} = [super{current_count} img(j,i)];
        % What are the current neighbours 
        [~, neighs] = find(Am(current_count,:) == 1);
        for ii = 1:(max(size(neighs)))
            neighbours{neighs(ii)} = [neighbours{neighs(ii)} img(j,i)];
        end
    end
end

%%%%
% TEST THROUGHOUT IMAGE
%%%%
for i = 1:h %swapped l and h around - NOTE!
    for j = 1:l
        cut = img(i,j);
        
        % Get Bands
        current_super = label_img(i,j);
        guardband = super{current_super};
        background = neighbours{current_super};
     

        % STATISTICAL CFAR (Weibull)
        pd = mle((background+1),'distribution','weibull');
        scale = pd(1,1);
        shape = pd(1,2);
        term = nthroot(-log(PFA),shape);
        TH = scale*term;

        % Check if CUT is target
        if cut >= TH 
            maskimg(i,j) = 1;
            detections = detections + 1;
        else 
            maskimg(i,j) = 0;
        end
        
    end
end
figure,imshow(maskimg,[]);






