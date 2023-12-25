% MATLAB R2018b
% Author: Yanting Yang yy3189@columbia.edu
% Modified Date: 2023-12-25

% Get all nii file pathes.
nii_4d_pathes = dir('./example_data/*/4d.nii.gz');

% Prepare pool for parallel execution.
% poolsize is recommended to be NI, which is 18 for the example data.
poolsize = 18;
parpool(poolsize);

% Process each nii file.
for i = 1:length(nii_4d_pathes)
    % Get nii path from index.
    nii_4d_path = fullfile(nii_4d_pathes(i).folder, nii_4d_pathes(i).name);
    
    % Print progress.
    fprintf('%d/%d\n', i, length(nii_4d_pathes));
    
    % Print nii path.
    disp(nii_4d_path);
    
    % Read nii file.
    nii_4d_data = double(niftiread(nii_4d_path));
    
    % Print nii data shape and type.
    whos nii_4d_data;
    
    % Get repetition, slice, length, and width from data shape.
    NR = size(nii_4d_data, 4);
    NI = size(nii_4d_data, 3);
    Ny = size(nii_4d_data, 1);
    Nx = size(nii_4d_data, 2);
    
    % Acquisition parameters
    r1 = 2.6;           % T1 relaxivity of Omniscan Gd-DTPA (mM-1 s-1)
    T10 = 0.9;          % T1 relaxation time of brain tissue before Gd-DTPA (s)
    TA = 33.6;          % Acquisition time of DCE-MRI sequence
    NormLimit = 0.1;    % Low norm -> better fit
    
    % algo: 1 for trust-region selective, 2 for levenberg-marquardt
    algo = 2;

    % AIF parameters
    A1 = 0.019556;
    m1 = -0.059213;
    A2 = 1.2419;
    m2 = 0.1689;
    
    % Image filter
    scont = zeros(size(nii_4d_data));
    fscont = zeros(size(nii_4d_data));
    for rep_idx = 1:NR
        for slc_idx = 1:NI
            scont(:,:,slc_idx,rep_idx) = nii_4d_data(:,:,slc_idx,rep_idx);
            fscont(:,:,slc_idx,rep_idx) = imgaussfilt(nii_4d_data(:,:,slc_idx,rep_idx),0.5);
        end
    end
    
    % Average over pre-contrast slices. For example data, pre-contrast
    % slices are slice 1 to slice 4.
    scont_pre = mean(scont(:,:,:,1:4),4);
    fscont_pre = mean(fscont(:,:,:,1:4),4);
    
    % Gd concentration
    Ce = zeros(size(nii_4d_data));
    for rep_idx = 1:NR
        Ce(:,:,:,rep_idx) = (fscont(:,:,:,rep_idx) - fscont_pre)./fscont_pre/r1/T10;
    end

    % Calculate over the whole image region.
    pos = [1 1 Nx-1 Ny-1];

    % Parallel execution
    Ktrans_map_full = zeros(size(scont_pre));
    parfor slc_idx = 1:NI
        % The funCalcCe will print at each line with format 
        % like (slc_idx, y_idx, 1).
        [Ktrans_map, Kep_map, t0_map, norm_map, count] = ...
            funCalcCe(Ce, NR, A1, A2, m1, m2, slc_idx, pos, algo, TA, NormLimit);
        Ktrans_map_full(:,:,slc_idx) = Ktrans_map(:,:,slc_idx);
    end
    
    % Save
    nii_info = niftiinfo(nii_4d_path);
    nii_info.BitsPerPixel = 64;
    nii_info.Datatype = 'double';
    nii_info.ImageSize = nii_info.ImageSize(1:3);
    nii_info.PixelDimensions = nii_info.PixelDimensions(1:3);
    nii_save_path = fullfile(nii_4d_pathes(i).folder, 'ktrans.nii');
    niftiwrite(Ktrans_map_full, nii_save_path, nii_info, 'Compressed', true);
end