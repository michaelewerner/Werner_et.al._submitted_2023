%% Run wavelet analysis and fit to an ANH model to analysis oscillatory signal
% Runs wavelet synchrosqueeze transform and fits ANH model to input data
% segmented and smoothed in python (see RingSegmentation and
% smoothing.ipynb) 
% Authors:
% Dylan Ray, Coleman Breen and Michael Werner
% BSD License
% December 2023

%% Save smoothed speed data to HDF5 file

c (folder name needs to be numerical) (Example: input = [1 2 3 4] where each numbered folder contains segmentation data for one ring )

n = size(input,1)

for i = 1:n
    k = input(i)
    filename = sprintf('input_file path/%d/SpeedDiffCS.xlsx',k) %%<< INSERT FOLDER PATH (example filename = sprintf('C:/data/%d/SpeedDiffCS.xlsx',k)
    speed = xlsread(filename,1)

    speedAdj = speed'

    %save adjusted speed data to HDF5 file
    filename2 = sprintf('/insert_folder_path/%d/speedAdj',k) %%<< INSERT FOLDER PATH in HDF5 file (filename2 = sprintf('/data/%d/speedAdj',k)
    s = size(speedAdj,1);
    h5create('HDF5_folder_path_and_name.hdf5',filename2, [s 72]) % << SPECIFY DESTINATION HDF5 FILE TO CREATE NEW FILE LOCATION (example: h5create('C:\data.hdf5',filename2, [s 72])
    h5write('HDF5_folder_path_and_name.hdf5',filename2, speedAdj) % << SPECIFY DESTINATION HDF5 FILE TO WRITE TO NEW FILE LOCATION (example : h5write('C:\data.hdf5',filename2, speedAdj)
end
%% Run wavelet synchrosqueeze transform analysis on all outputs from any given condition
input = [] %% <<< list input folders (folder name needs to be numerical)

Recon_by_Segcat = []     
pRidgeAllcat = []
IpsAllcat = []
AmpsAllcat = []
histANHcat = []
histANHadj2cat = []
vintervalANHcat = []
vintervalANHadj2cat = []
histANHAdjcat = []
histANHadj2cat = []

% Default parameters for wavelet computation
freqBins = 16;
nRidges = 25;
nPixels = 16;
penalty = 2;

n = size(input)

for i = 1:n
    k = input(i)
    filename = sprintf('/input_file_path/%d/speedAdj',k) %<<< specify input file path (example: h5create (filename = sprintf(('/data/%d/speedAdj',k))
    speedAdj = h5read('specify_filepath_and_HDF5_file_name.hdf5', filename) %<<< specify input file path (example: speedAdj = h5read('C:\data.hdf5', filename)
    
    signals = speedAdj
    signals = diff(signals)
    s= size(signals,1)
    signals = padarray(signals,s,0,'both')
    
    ts = seconds(2.7); % frames recorded every 2.7 secs
    signals = double(signals)
   
    % Compute the wsst for each element in this sample
    [sst_single, ~] = wsst(signals(:,1), ts, 'amor', 'VoicesPerOctave', 32); % compute on one element
    nrows = size(sst_single, 1);
    ncols = size(sst_single, 2);


    sst_by_time_real = NaN(nrows, ncols, 72);
    sst_by_time_imag = NaN(nrows, ncols, 72);
    sst_by_time = NaN(nrows, ncols, 72);

    frequency_by_time = NaN(nrows, 1, 72);

    pridge_by_time = NaN(size(signals, 1), nRidges, 72);
    rec_modes = NaN(size(signals, 1), nRidges, 72);
    iridge_by_time = NaN(size(signals, 1), nRidges, 72);
    filtered_signals = NaN(size(signals));
    input_signal = []
    
    for u = 1:72
          one_signal = signals(:,u);


          [sst_temp, frequency] = wsst(one_signal, ts, 'amor', 'VoicesPerOctave', 32);


          input_signal = [input_signal one_signal]


          period = 1./frequency;

          period_matrix = repmat(period', 1, ncols);

    %       Optional upper and lower cutoffs
    %       sst_temp(period_matrix > 120) = 0;
    %       sst_temp(period_matrix < 15) = 0;

          filtered_signals(:,u) = iwsst(sst_temp);
 

          [fridge,iridge] = wsstridge(sst_temp,penalty, frequency, 'NumRidges', nRidges,'NumFrequencyBins', nPixels);
          pridge = 1./fridge;   

          xrec = iwsst(sst_temp, iridge,'NumFrequencyBins', freqBins); %,'NumFrequencyBins', freqBins
      

          % Keep absolute value (eliminates complex part and negatives)
          sst_by_time_real(:,:,u) = real(sst_temp);
          sst_by_time_imag(:,:,u) = imag(sst_temp);

          sst_by_time(:,:,u) = abs(sst_temp);

      
          % Keep a matrix of all of the frequency levels
          frequency_by_time(:,:,u) = frequency;
         

          % Period ridge
           pridge_by_time(:,:,u) = pridge;


           % Index ridge
           iridge_by_time(:,:,u) = iridge;


           % Reconstructed modes
           rec_modes(:,:,u) = xrec;

     end


%
    % function fit_ANH(ifile, ofile, max_num_modes)
    %FIT_ANH fits the adaptive non-harmonic model to data.

    % Daubechies, Ingrid & Lu, Jianfeng & Wu, Hau-Tieng. (2011).
    % Synchrosqueezed wavelet transforms: An empirical mode decomposition-like tool.
    % Applied and Computational Harmonic Analysis. 30. 243-261. 10.1016/j.acha.2010.08.002.

    % Initialize root mean squared error (cost function)
       RMSE = zeros(72, 1);
       max_num_modes = 5

    IMFs_by_Seg = []
    amps_by_Seg = []
    ips_by_Seg = []
    Recon_by_Seg = []
    ReconAmps_by_Seg = []
    ReconIps_by_Seg = []
    mIndex_by_Seg = []
 
    for j = 1:72
           signal = signals(:, j);
           rec = rec_modes(:, :, j);
           ridge = pridge_by_time(:, :, j);

           [RMSE(j), mindex] = find_best_IMFs(signal, rec, max_num_modes);

           IMFs = rec(:, mindex); % intrinsic mode functions
           amps = abs(hilbert(IMFs)); % amplitudes of corresponding modes
           ips = ridge(:, mindex); % instantaneous period 
           Recon = sum(IMFs,2)
           ReconAmps = sum(amps,2)
           ReconIps = sum(ips,2)


           s1 = size (IMFs,1)
           sz = size(IMFs,2)
           s2 = max_num_modes - sz
           sIdx = size(mindex,2)
           if sz == max_num_modes
           elseif sz < max_num_modes 
               IMFs = [IMFs nan(s1,s2)]
               amps = [amps nan(s1,s2)]
               ips = [ips nan(s1,s2)]
               mindex = [mindex nan(1,s2)]
           end
           
           mIndex_by_Seg(:,:,j) = mindex
           IMFs_by_Seg(:,:,j) = IMFs;
           amps_by_Seg(:,:,j) = amps;
           ips_by_Seg(:,:,j) = ips;
           Recon_by_Seg(:,j) = Recon;
           ReconAmps_by_Seg(:,j) = ReconAmps;
           ReconIps_by_Seg(:,j) = ReconIps;

    end
    
   
    s=size(signals,1)
    
    %Concatenate all instateneous periods and amplitudes
    IpsAll = []
    for i = 1:5
        cat = reshape(ips_by_Seg(s:end-s,i,:),[],1)
        IpsAll = [IpsAll; cat]
    end

    AmpsAll = []
    for i = 1:5
        cat = reshape(amps_by_Seg(s:end-s,i,:),[],1)
        AmpsAll = [AmpsAll; cat]
    end
    
    % save all padded inputs to HDF5
    s=size(signals,1)
    filename = sprintf('/insert_folder_path/%d/padded_input',k) %%<< INSERT FOLDER PATH in HDF5 file (filename = sprintf('/data/%d/padded_input',k)
    h5create('HDF5_folder_path_and_name.hdf5', filename, [s 72]) %% << SPECIFY DESTINATION HDF5 FILE TO CREATE NEW FILE LOCATION (example: h5create('C:\data.hdf5',filename, [s 72])
    h5write('HDF5_folder_path_and_name.hdf5', filename, signals) % << SPECIFY DESTINATION HDF5 FILE TO WRITE TO NEW FILE LOCATION (example : h5write('C:\data.hdf5',filename, padded_input)
    
    %Save concatenated instanteneous periods
    s=size(IpsAll,1)
    filename = sprintf('/insert_folder_path/%d/IpsAll',k) %%<< INSERT FOLDER PATH in HDF5 file (filename = sprintf('/data/%d/IpsAll',k)
    h5create('HDF5_folder_path_and_name.hdf5', filename, [s]) %% <<SPECIFY DESTINATION HDF5 FILE TO CREATE NEW FILE LOCATION (example: h5create('C:\data.hdf5',filename, [s])
    h5write('HDF5_folder_path_and_name.hdf5', filename, IpsAll) % << SPECIFY DESTINATION HDF5 FILE TO WRITE TO NEW FILE LOCATION (example : h5write('C:\data.hdf5',filename, IpsAll)
    
    %Save concatenated Amplitudes
    s=size(AmpsAll,1)
    filename = sprintf('/insert_folder_path/%d/AmpsAll',k) %%<< INSERT FOLDER PATH in HDF5 file (filename = sprintf('/data/%d/AmpsAll',k)
    h5create('HDF5_folder_path_and_name.hdf5', filename, [s]) %% <<SPECIFY DESTINATION HDF5 FILE TO CREATE NEW FILE LOCATION (example: h5create('C:\data.hdf5',filename2, [s])
    h5write('HDF5_folder_path_and_name.hdf5', filename, AmpsAll) % << SPECIFY DESTINATION HDF5 FILE TO WRITE TO NEW FILE LOCATION (example : h5write('C:\data.hdf5',filename2, AmpsAll)
 
end

%% Plot Weighted histograms aggregating all mode frequencies from the ANH fit with the lowest RMSE adjusting the amplitudes by their corresponding angular frequencies 

input = []

IpsAll = []
AmpsAll = []

n = size(input,1)

%% Load and Concatenate all instantenous periods and amplitudes for all data from a single condition
for i = 1:n
    k = input(i)
    filename = sprintf('/insert_folder_path/%d/IpsAll',k) %<<< specify input file path (example: h5create (filename = sprintf(('/data/%d/IpsAll',k))
    ips = h5read('HDF5_folder_path_and_name.hdf5', filename) %<<< specify input file path (example: ips = h5read('C:\data.hdf5', filename)
    
    filename = sprintf('/insert_folder_path/%d/AmpsAll',k) %<<< specify input file path (example: h5create (filename = sprintf(('/data/%d/AmpsAll',k))
    amps = h5read('HDF5_folder_path_and_name.hdf5', filename) %<<< specify input file path (example: amps = h5read('C:\data.hdf5', filename)
    
    IpsAll = [IpsAll; ips]
    AmpsAll = [AmpsAll; amps]
end 

IpsAlln = IpsAll(~isnan(IpsAll))
AmpsAlln = AmpsAll(~isnan(AmpsAll))
s3=round(max(IpsAlln))

%adjust amplitudes by their corresponding angular frequencies and bin 
AmpsAllAdj = AmpsAlln./(2*pi*(1./IpsAlln))
[ampsANH, intANH] = histwc(IpsAll, AmpsAll, s3)
[histANHadj2, vintervalANHadj2] = histwc(IpsAlln, AmpsAllAdj, s3)

% plot weighted histogramm
figure
bar(vintervalANHadj2,histANHadj2,'histc')
xlim([0 200])
xlabel('period (seconds)')
ylabel('relative distribution')


%%%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [min_rmse, mindex] = find_best_IMFs(signal, rec, max_nm, min_nm)
%signal: zero centered signal (inward velocity of one of 72 elements
%towards center)
%rec: reconstuction components (AM-FM basis functions). These are the
%pieces we use to rebuild the signal. In the wavelet computation, we picked
%the 25 best recs, and now we'll perform best subset selection to get 3,4,5
%that represent the signal well
%max_nm: max_number of modes
%min_nm: min number of modes (set this = to max_nm if you want exactly that
%many of modes)

    % Vestige from experimenting. Ignore but NTS: don't delete
    if nargin == 3
        min_nm = 1;
    end

    num_modes = size(rec, 2);
    num_combos = 0;
    
    % Number of combinations w/ 1, 2, ..., modes out of all possible (prob
    % 25)
    for m = min_nm:max_nm
        num_combos = num_combos + nchoosek(num_modes, m);
    end
    
    % initalize matrix
    all_combos = nan(num_combos, max_nm);
    ix = 1;
    
    for m = min_nm:max_nm
       temp = combnk(1:num_modes, m);
       
       nrow = length(temp);
       ncol = size(temp, 2);
       
       all_combos(ix:ix+nrow-1, 1:ncol) = temp;
       ix = ix + nrow;
    end
    
    % Relative Root Mean Squared Error
    RMSE = zeros(num_combos, 1);
    
    for combo = 1:num_combos
        
        combo_ix = all_combos(combo, :); % pick one set of indicies (e.g. [4, 6, 15, NaN, NaN])
        combo_ix(isnan(combo_ix)) = []; % lop off the nans
        
        temp_modes = rec(:, combo_ix); % grab corresponding modes
        recon = sum(temp_modes, 2); % reconstruct signal with specified modes
        
        % Relative Root Mean Squared Error
        RMSE(combo) = sqrt( mean((signal - recon).^2 )) / (max(signal) - min(signal));
    end
    
    % Find the best combination after having computed all
    [min_rmse, combo_mindex] = min(RMSE);
    mindex = all_combos(combo_mindex, :);
    
    mindex(isnan(mindex)) = [];
end

% HISTWC  Weighted histogram count given number of bins
%
% This function generates a vector of cumulative weights for data
% histogram. Equal number of bins will be considered using minimum and 
% maximum values of the data. Weights will be summed in the given bin.
%
% Usage: [histw, vinterval] = histwc(vv, ww, nbins)
%
% Arguments:
%       vv    - values as a vector
%       ww    - weights as a vector
%       nbins - number of bins
%
% Returns:
%       histw     - weighted histogram
%       vinterval - intervals used
%       
%
%
% See also: HISTC, HISTWCV
% Author:
% mehmet.suzen physics org
% BSD License
% July 2013

function [histw, vinterval] = histwc(vv, ww, nbins)
  minV  = min(vv);
  maxV  = max(vv);
  delta = (maxV-minV)/nbins;
  vinterval = linspace(minV, maxV, nbins)-delta/2.0;
  histw = zeros(nbins, 1);
  for i=1:length(vv)
    ind = find(vinterval < vv(i), 1, 'last' );
    if ~isempty(ind)
      histw(ind) = histw(ind) + ww(i);
    end
  end
end