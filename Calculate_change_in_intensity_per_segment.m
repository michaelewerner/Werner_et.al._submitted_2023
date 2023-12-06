%% extract fluorescent intensity data from ring segments
% Calculated average intensity data for each ring segment from segment coordinates determined using Python script to segment ring (see RingSegmentation and
% smoothing.ipynb) by averaging the intensity of each poligine segment and
% the 8 surounding pixels.
% Author:
% Michael Werner
% BSD License
% December 2023

%Import intensity data from input Tif image.
tiff_info = imfinfo('Insert input file path and filename.tif'); % <<< INSERT INPUT FILEPATH AND FILENAME
tiff_stack = imread('Insert input file path and filename', 1) ;  % <<< INSERT INPUT FILEPATH AND FILENAME
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread('Insert input file path and filename', ii); % <<< INSERT INPUT FILEPATH AND FILENAME
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end

%Import smoothed coordinates of ringe segments 
XnewCS=xlsread('Insert Filepath and Filename.xlsx',1) % <<< INSERT INPUT FILEPATH AND FILENAME of XYs.xlsx file created using RingSegmentation and
% smoothing.ipynb
YnewCS=xlsread('Insert Filepath and Filename.xlsx',2) % <<< INSERT INPUT FILEPATH AND FILENAME of XYs.xlsx file created using RingSegmentation and
% smoothing.ipynb

%calculate avereage intensity for each segment for each timepoint
s = size (tiff_stack,3)
IntAll = []
IntMean = []
for t = 1:s
    IntA = []
    IntM = []
    S = squeeze(tiff_stack(:,:,t)')
    for c = 1:72
        x = round(XnewCS(c,t)) 
        y = round(YnewCS(c,t))
        
        x1 = round(XnewCS(c,t)) + 1
        y1 = round(YnewCS(c,t)) + 1
        
        x2 = round(XnewCS(c,t)) - 1
        y2 = round(YnewCS(c,t)) - 1
       
        int1 = S(x,y)
        int2 = S(x,y1)
        int3 = S(x,y2)
        int4 = S(x1,y)
        int5 = S(x1,y1)
        int6 = S(x1,y2)
        int7 = S(x2,y)
        int8 = S(x2,y1)
        int9 = S(x2,y2)
        
        intA = [int1 int2 int3 int4 int5 int6 int7 int8 int9]
        intM = nanmean(intA)

        IntA = [IntA; intA]
        IntM = [IntM intM]
    end
    IntA = double (IntA)
    IntAll (:,:,t) = IntA
    IntMean = [IntMean; IntM]
end

% Calculate smoothed intensity
IntMeanMov = movmean(IntMean,5,1)
IntMeanMov2 = movmean(IntMeanMov,5,2)

% Calculate smoothed change in intensity
IntMeanDiff = diff(IntMean)
IntMeanDiffmov = movmean(IntMeanDiff,5,1)
IntMeanDiffmov2 = movmean(IntMeanDiffmov,5,2)