function varargout = Permeability_estimation(varargin)
% PERMEABILITY_ESTIMATION MATLAB code for Permeability_estimation.fig
%      PERMEABILITY_ESTIMATION, by itself, creates a new PERMEABILITY_ESTIMATION or raises the existing
%      singleton*.
%
%      H = PERMEABILITY_ESTIMATION returns the handle to a new PERMEABILITY_ESTIMATION or the handle to
%      the existing singleton*.
%
%      PERMEABILITY_ESTIMATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PERMEABILITY_ESTIMATION.M with the given input arguments.
%
%      PERMEABILITY_ESTIMATION('Property','Value',...) creates a new PERMEABILITY_ESTIMATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Permeability_estimation_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Permeability_estimation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Permeability_estimation

% Last Modified by GUIDE v2.5 06-Jul-2023 07:16:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Permeability_estimation_OpeningFcn, ...
                   'gui_OutputFcn',  @Permeability_estimation_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end 

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Permeability_estimation is made visible.
function Permeability_estimation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Permeability_estimation (see VARARGIN)

% Choose default command line output for Permeability_estimation
handles.output = hObject;

% Initialize variables

handles.DCE_acquisitions = 64; % 84 Numbe of DCE acquisitions
handles.DCE_slices = 18; % 18 Number of DCE slices
handles.acquisition_time = 33.6; %22 38; % Acquisition time (TA) (seconds)
handles.repetition_time = 0.2; %0.2 0.23; % Repetition time (TR)  (seconds)

handles.Gd_relaxivity = 2.6; % T1 relaxivity of Omniscan Gd-DTPA (mM-1 s-1)
handles.T1_brain_relaxation_time = 0.9; % T1 relaxation time of brain tissue (seconds)
handles.T1_blood_relaxation_time = 1.5; % T1 relaxation time of blood (seconds)
handles.normalization_limit = 0.1;

handles.current_MRI_slice = 1;
handles.current_MRI_acquisition = 1;
handles.current_Ktrans_slice = 18;
handles.current_Ktrans_acquisition = 64;
handles.Ktrans_acquisitions = 64;
handles.Ktrans_slices = 18;

handles.ROI_position = [10 10 140 140];
handles.selection_preference = 1; % 1: Analyze current MRI slice and acquisition, 2: analyze selected range of MRI slices and acqusitions

handles.max_Ktrans = 0.1;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Permeability_estimation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Permeability_estimation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function current_MRI_slice_slider_Callback(hObject, eventdata, handles)
% hObject    handle to current_MRI_slice_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.current_MRI_slice = uint8(get(hObject, 'Value'));

current_slice = handles.current_MRI_slice;
current_acquisition = handles.current_MRI_acquisition;

axes(handles.axes1);
imagesc(handles.DCE_MRI_frames{current_acquisition,current_slice}); colormap(gray);
%imagesc(handles.DCE_MRI_frames{current_acquisition,current_slice},[0 256]); colormap(gray);
xlabel('x (pixels)'); ylabel('y (pixels)'); C=colorbar; ylim(C,[0 256]);

delete(handles.slice_text);
slice_number = ['Slice: ',num2str(current_slice),'/',num2str(handles.DCE_slices)];
handles.slice_text = annotation('textbox',[0.34 0.85 0.1 0.1],'String',slice_number,'FontSize',14,'Color','w','LineStyle','none');

guidata(hObject, handles); % Update the GUI data structure


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.current_Ktrans_slice = uint8(get(hObject, 'Value'));

current_Ktrans_slice=handles.current_Ktrans_slice;
current_Ktrans_acquisition=handles.current_Ktrans_acquisition;

Ktrans_map = handles.Ktrans_map{current_Ktrans_acquisition,current_Ktrans_slice};
acquisition_time = handles.acquisition_time;
max_Ktrans = handles.max_Ktrans;

Ktrans_img = double(Ktrans_map(:,:)/acquisition_time*60); % min-1
 
axes(handles.axes2); colormap(handles.axes2,jet)
imagesc(Ktrans_img,[0 handles.max_Ktrans]); C=colorbar; 
ylim(C,[0 max_Ktrans]);
xlabel('x (pixels)'); ylabel('y (pixels)'); ylabel(C,'Ktrans (min-1)');

delete(handles.slice_text_2);
slice_number_2 = ['Slice: ',num2str(handles.current_Ktrans_slice),'/',num2str(handles.DCE_slices)];
handles.slice_text_2 = annotation('textbox',[0.68 0.85 0.1 0.1],'String',slice_number_2,'FontSize',14,'Color','w','LineStyle','none');

guidata(hObject, handles); % Update the GUI data structure


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure_name = ['MRI_slice_',num2str(handles.current_MRI_slice),'_acquisition_',num2str(handles.current_MRI_acquisition)];

fig = figure; 
copyobj(handles.axes1,fig); colormap(gray); daspect([1 1 1]);
hgsave(fig,[figure_name,'.fig']);
saveas(fig,[figure_name,'.png']);

delete(fig);

guidata(hObject, handles); % Update the GUI data structure


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure_name = ['Ktrans_slice_',num2str(handles.current_Ktrans_slice),'_acquisition_',num2str(handles.current_Ktrans_acquisition)];

fig = figure; 
copyobj(handles.axes2,fig); colormap(jet); daspect([1 1 1]);
hgsave(fig,[figure_name,'.fig']);
saveas(fig,[figure_name,'.png']);

delete(fig);

guidata(hObject, handles); % Update the GUI data structure



% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[handles.x_selected,handles.y_selected]=ginput(1);
% disp([handles.x_selected,handles.y_selected]);

axes(handles.axes3);
x_point = int8(round(handles.x_selected));
y_point = int8(round(handles.y_selected));
disp([x_point,y_point]);

Ktrans_per_point=zeros(1,handles.DCE_acquisitions);
current_Ktrans_slice = handles.current_Ktrans_slice;
selected_initial_acquisition = handles.selected_initial_acquisition;
selected_final_acquisition = handles.selected_final_acquisition;

for i =  selected_initial_acquisition:selected_final_acquisition
    
    Ktrans_acquisition = handles.Ktrans_map{i,current_Ktrans_slice};
    Ktrans_per_point(1,i) = Ktrans_acquisition(y_point,x_point);
    
end

plot((1:length(Ktrans_per_point))*handles.acquisition_time/60,Ktrans_per_point/handles.acquisition_time*60,'o','Color','r');
ylim([0 handles.max_Ktrans]);
xlabel('Time (min)'); ylabel('Ktrans (min-1)');

set(handles.Ktrans_point_coordinates,'String',['x: ',num2str(x_point),' pixels',sprintf('\n'),'y: ',num2str(y_point),' pixels']);

handles.Ktrans_per_point = Ktrans_per_point;

guidata(hObject, handles); % Update the GUI data structure



% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure_name = ['Ktrans_over_time_slice_',num2str(handles.current_Ktrans_slice),'_x_',num2str(round(handles.x_selected)),'_y_',num2str(round(handles.y_selected))];

fig = figure; 
copyobj(handles.axes3,fig);
hgsave(fig,[figure_name,'.fig']);
saveas(fig,[figure_name,'.png']);

delete(fig);

guidata(hObject, handles); % Update the GUI data structure



% --- Executes on button press in load_data.
function load_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% [MRI_frames] = uigetfile({'*.dcm','DICOM files (*.dcm)';'*.jpeg;*.png;*.bmp','Image files (*.jpeg,*.png,*.bmp)'},'Open DCE-MRI file','Multiselect','on');
% [MRI_frames] = imgetfile('Multiselect','on');
addpath(genpath('/home/sail/matlab/matlab_toolbox/'))
addpath(genpath('/home/sail/matlab/matlab_toolbox/NIfTI_20140122'))
addpath(genpath('matlab_toolbox'))

[MRI_frames, pathname] = uigetfile({'*.nii.gz'});
MRI_4d_str = load_nii(fullfile(pathname, MRI_frames));
MRI_4d_volume = double(MRI_4d_str.img);
DCE_MRI_frames = cell(handles.DCE_acquisitions,handles.DCE_slices);

% Reordering % 

for i = 1:handles.DCE_acquisitions
    
    for j = 1:handles.DCE_slices
    
        selected_frame = (i-1)*handles.DCE_slices+j;
        DCE_MRI_frames{i,j} = MRI_4d_volume(:,:,j,i);
        %DCE_MRI_frames{i,j} = imread(MRI_frames{selected_frame});
    
    end
    
end

current_slice = 1;
current_acquisition = 1;

axes(handles.axes1);
imagesc(DCE_MRI_frames{current_acquisition,current_slice}); colormap(gray);
%imagesc(DCE_MRI_frames{current_acquisition,current_slice},[0 256]); colormap(gray);
xlabel('x (pixels)'); ylabel('y (pixels)'); C=colorbar; ylim(C,[0 256]);

slice_number = ['Slice: ',num2str(current_slice),'/',num2str(handles.DCE_slices)];
handles.slice_text = annotation('textbox',[0.34 0.85 0.1 0.1],'String',slice_number,'FontSize',14,'Color','w','LineStyle','none');
acquisition_number = ['Acquisition: ',num2str(current_acquisition),'/',num2str(handles.DCE_acquisitions)];
handles.acquisition_text = annotation('textbox',[0.34 0.83 0.1 0.1],'String',acquisition_number,'FontSize',14,'Color','w','LineStyle','none');

handles.DCE_MRI_frames = DCE_MRI_frames;
handles.pathname = pathname;
guidata(hObject, handles); % Update the GUI data structure



% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes2);
% running_text = annotation('textbox',[0.68 0.85 0.1 0.1],'String','Running...','FontSize',14,'Color','k','LineStyle','none');
delete(handles.estimated_time_text);

NR = handles.DCE_acquisitions;
r1 = handles.Gd_relaxivity;                 %T1 relaxivity of Omniscan Gd-DTPA (mM-1 s-1)
T10 = handles.T1_brain_relaxation_time;     %T1 relaxation time of brain tissue before Gd-DTPA (s)
T10art = handles.T1_blood_relaxation_time;  %T1 relaxation time of blood before Gd-DTPA (s)
TR = handles.repetition_time;               %Repetition time of DCE-MRI sequence
TA = handles.acquisition_time;              %Acquisition time of DCE-MRI sequence
Hct=0.45;                                   %Hematocrit
FilterMatrix = ones(3,3) / 9;               %Smoothing pixel data with neighbors
NormLimit = 0.1;                            % Low norm -> better fit

% Algorithms used for model fitting: 
% 1 for trust-region reflective
% 2 for levenberg-marquardt
algoAIF = 2;    % fit for AIF
algoCe = 2;     % fit for pixel-by-pixel estimation
algoCeROI = 2;  % fit for quantitative measurement

% AIF parameters
A1 = 0.019556;
m1 = -0.059213;
A2 = 1.2419;
m2 = 0.1689;

% Non-linear least squares method for Ktrans with AIF model fitting.

%  Pre-allocate matrices % 

Ktrans_map_full = cell(handles.DCE_acquisitions,handles.DCE_slices);
Kep_map_full = cell(handles.DCE_acquisitions,handles.DCE_slices);
t0_map_full = cell(handles.DCE_acquisitions,handles.DCE_slices);
norm_map_full = cell(handles.DCE_acquisitions,handles.DCE_slices);
count_full = cell(handles.DCE_acquisitions,handles.DCE_slices);

position = handles.ROI_position;
current_MRI_acquisition = double(handles.current_MRI_acquisition);
current_MRI_slice = handles.current_MRI_slice;

% % Calculate Ktrans map only for current MRI slice and acquisition % %

if handles.selection_preference == 1 
    acquisition = current_MRI_acquisition; % Must be double!
    %sel_slice = current_MRI_slice;
    
    % Filter images %
    
    for i=1:acquisition
        for j=1:handles.DCE_slices        
            scont(:,:,j,i)=handles.DCE_MRI_frames{i,j};
            %fscont(:,:,j,i)=imfilter(handles.DCE_MRI_frames{i,j},FilterMatrix); % Apply filter
            fscont(:,:,j,i)=imgaussfilt(handles.DCE_MRI_frames{i,j},0.5);
        end
    end

    % Pre-contrast %

    for j=1:handles.DCE_slices
        scont_pre(:,:,j)=handles.DCE_MRI_frames{1,j};
        %fscont_pre(:,:,j)=imfilter(handles.DCE_MRI_frames{1,j},FilterMatrix); % Apply filter
        fscont_pre(:,:,j)=imgaussfilt(handles.DCE_MRI_frames{1,j},0.5); % Apply filter
    end

    fscont = double(fscont);        % smoothed tissue pixels
    scont = double(scont);          % original tissue pixels
    scont_pre = double(scont_pre);
    fscont_pre = double(fscont_pre);

    [Ny, Nx, Ns, NI]=size(fscont);

    Ce = zeros(Ny,Nx,Ns,acquisition);

        for nI = 1:acquisition
            Ce(:,:,:,nI) = (fscont(:,:,:,nI) - ...
                fscont_pre(:,:,:))./(fscont_pre(:,:,:))/r1/T10;
        end

    NormLimit = 0.1; % Low norm -> better fit 
    Ktrans_res = zeros(size(scont_pre));
    Kep_res = zeros(size(scont_pre));
    
    %%0506
    for acquisition = [64, 84]
    for sel_slice=1:handles.DCE_slices
        [Ktrans_map, Kep_map, t0_map, norm_map, count] = funCalcCe(Ce, acquisition, ...
         A1, A2, m1, m2, sel_slice, position, algoCe, TA, NormLimit);

        Ktrans_map = squeeze(Ktrans_map(:,:,sel_slice));
        Kep_map = squeeze(Kep_map(:,:,sel_slice));

        slice = sel_slice;

        Ktrans_map_full{acquisition,slice} = Ktrans_map;
        Kep_map_full{acquisition,slice} = Kep_map;
        t0_map_full{acquisition,slice} = t0_map;
        norm_map_full{acquisition,slice} = norm_map;
        count_full{acquisition,slice} = count;
        Ktrans_res(:,:,slice) = Ktrans_map;
        Kep_res(:,:,slice) = Kep_map;
    end
%     set(handles.slider1,'Visible','Off');
%     set(handles.slider2,'Visible','Off');

    handles.current_Ktrans_acquisition = handles.current_MRI_acquisition;
    handles.current_Ktrans_slice = handles.current_MRI_slice;
    
       save(fullfile(handles.pathname, ['Ktrans_map_', num2str(acquisition), '.mat']), 'Ktrans_res', '-mat');
       save(fullfile(handles.pathname, ['Kep_map_', num2str(acquisition), '.mat']), 'Kep_res', '-mat');    
    end
    save(fullfile(handles.pathname, 'Ktrans_map_full.mat'), 'Ktrans_map_full', '-mat');

    % % Calculate Ktrans map for selected MRI slices and acquisitions % %
     
elseif handles.selection_preference == 2 % Calculate Ktrans maps for all selected slices/acquisitions

    slice = 0;
    h = waitbar(0,['Please wait... Completion: ',num2str(100*slice/handles.DCE_slices,2),'%']);
    
    selected_initial_slice = double(handles.selected_initial_slice);
    selected_final_slice = double(handles.selected_final_slice);
    selected_reference_acquisition = double(handles.selected_reference_acquisition);
    selected_initial_acquisition = double(handles.selected_initial_acquisition);
    selected_final_acquisition = double(handles.selected_final_acquisition);
   
    s=1;   
    
    for i=1:selected_final_acquisition
        for j=1:handles.DCE_slices        
            scont(:,:,j,i)=handles.DCE_MRI_frames{i,j};
            %fscont(:,:,j,i)=imfilter(handles.DCE_MRI_frames{i,j},FilterMatrix);
            fscont(:,:,j,i)=imgaussfilt(handles.DCE_MRI_frames{i,j},0.5);
        end
    end
    
    % Pre-contrast %
    scont_pre=mean(scont(:,:,:,1:selected_reference_acquisition),4);
    fscont_pre=mean(fscont(:,:,:,1:selected_reference_acquisition),4);
    %for j=1:handles.DCE_slices
    %    scont_pre(:,:,j)=handles.DCE_MRI_frames{1,j};
    %    fscont_pre(:,:,j)=imfilter(handles.DCE_MRI_frames{1,j},FilterMatrix); % Apply filter
    %end
    
    fscont = double(fscont);        % smoothed tissue pixels
	scont = double(scont);          % original tissue pixels
	scont_pre = double(scont_pre);
	fscont_pre = double(fscont_pre);
            
    Ktrans_res = zeros(size(scont_pre));
    Kep_res = zeros(size(scont_pre));
    s = 1;
    for acquisition = selected_final_acquisition
%0505    for acquisition = selected_initial_acquisition:selected_final_acquisition
       for slice = selected_initial_slice:selected_final_slice
            [Ny, Nx, Ns, NI]=size(fscont);
            Ce = zeros(Ny,Nx,Ns,acquisition);
            sel_slice = slice;

            for nI = 1:acquisition
                Ce(:,:,:,nI) = (fscont(:,:,:,nI) - ...
                    fscont_pre(:,:,:))./(fscont_pre(:,:,:))/r1/T10;
            end

            NormLimit = 0.1; % Low norm -> better fit 
        
            [Ktrans_map, Kep_map, t0_map, norm_map, count] = funCalcCe(Ce, acquisition, ...
            A1, A2, m1, m2, sel_slice, position, algoCe, TA, NormLimit);
        
            Ktrans_map = squeeze(Ktrans_map(:,:,sel_slice));
            Kep_map = squeeze(Kep_map(:,:,sel_slice));
 
            Ktrans_map_full{acquisition,slice} = Ktrans_map;
            Kep_map_full{acquisition,slice} = Kep_map;
            Ktrans_res(:,:,slice) = Ktrans_map;
            Kep_res(:,:,slice) = Kep_map;
            t0_map_full{acquisition,slice} = t0_map;
            norm_map_full{acquisition,slice} = norm_map;
            count_full{acquisition,slice} = count;

       end
       save(fullfile(handles.pathname, ['Ktrans_map_', num2str(acquisition), '.mat']), 'Ktrans_res', '-mat');
       save(fullfile(handles.pathname, ['Kep_map_', num2str(acquisition), '.mat']), 'Kep_res', '-mat');    
    end

    delete(h);
    h = waitbar(s/(selected_final_acquisition-selected_initial_acquisition+1),['Please wait... Completion: ',num2str(100*s/(selected_final_acquisition-selected_initial_acquisition+1),2),'%']);
    set(h,'windowstyle','modal');
    s=s+1;
    delete(h);
%     set(handles.slider1,'Visible','On');
%     set(handles.slider2,'Visible','On');
    
end
 

 handles.fitting_results = {Ktrans_map_full,Kep_map_full,t0_map_full,norm_map_full}; 

 current_MRI_acquisition = handles.current_MRI_acquisition;
 current_MRI_slice = handles.current_MRI_slice;
 handles.current_Ktrans_slice = current_MRI_slice;
 %('Ktrans_map_full.mat', 'Ktrans_map_full', '-mat');

 Ktrans_plot = Ktrans_map_full{current_MRI_acquisition,current_MRI_slice};
 Ktrans_img = double(Ktrans_map(:,:)/TA*60); % min-1
 
 axes(handles.axes2); colormap(handles.axes2,jet);
 imagesc(Ktrans_img,[0 handles.max_Ktrans]); C=colorbar; 
 caxis([0 handles.max_Ktrans])
 xlabel('x (pixels)'); ylabel('y (pixels)'); ylabel(C,'Ktrans (min-1)');
 
 slice_number_2 = ['Slice: ',num2str(handles.current_MRI_slice),'/',num2str(handles.DCE_slices)];
 handles.slice_text_2 = annotation('textbox',[0.68 0.85 0.1 0.1],'String',slice_number_2,'FontSize',14,'Color','w','LineStyle','none');
 acquisition_number_2 = ['Acquisition: ',num2str(handles.current_MRI_acquisition),'/',num2str(handles.DCE_acquisitions)];
 handles.acquisition_text_2 = annotation('textbox',[0.68 0.83 0.1 0.1],'String',acquisition_number_2,'FontSize',14,'Color','w','LineStyle','none');
 
 handles.Ktrans_map = Ktrans_map_full;
 size(Ktrans_map_full)
 handles.Kep_map = Kep_map_full;
 handles.t0_map = t0_map_full;
 handles.norm_map = norm_map_full;
 handles.count = count_full;
 
 guidata(hObject, handles); % Update the GUI data structure


% --- Executes on button press in select_ROI.
function select_ROI_Callback(hObject, eventdata, handles)
% hObject    handle to select_ROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.estimated_time_text = annotation('textbox',[0.34 0.40 0.1 0.1],'String','Select ROI','FontSize',14,'Color','w','LineStyle','none');

%    ROI_map = [40 40 80 80];
    ROI_map = [1 1 200 211];
    h2 = imrect(gca,ROI_map);
    handles.ROI_position = wait(h2);
    disp(handles.ROI_position);
    
    number_of_pixels = round(handles.ROI_position(3))*round(handles.ROI_position(4)); % Number of pixels within the selected ROI
    estimated_time = number_of_pixels/45; % Estimated time in seconds (calculation speed: 45 pixels/sec)
    
    delete(handles.estimated_time_text);
    time_for_Ktrans_model = ['Estimated time: ',num2str(estimated_time,2),'s or ',num2str(estimated_time/60,2),'min'];
    handles.estimated_time_text = annotation('textbox',[0.34 0.40 0.1 0.1],'String',time_for_Ktrans_model,'FontSize',14,'Color','w','LineStyle','none');

guidata(hObject, handles); % Update the GUI data structure


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

handles.normalization_limit = str2double(get(hObject,'String'));

guidata(hObject, handles); % Update the GUI data structure


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double

handles.max_Ktrans = str2double(get(hObject,'String'));

current_Ktrans_slice=handles.current_Ktrans_slice;
current_Ktrans_acquisition=handles.current_Ktrans_acquisition;

Ktrans_map = handles.Ktrans_map{current_Ktrans_acquisition,current_Ktrans_slice};
acquisition_time = handles.acquisition_time;
max_Ktrans = handles.max_Ktrans;

Ktrans_img = Ktrans_map(:,:)/acquisition_time*60; % min-1
 
axes(handles.axes2); colormap(handles.axes2,jet);
imagesc(Ktrans_img,[0 max_Ktrans]); C=colorbar; 
ylim(C,[0 max_Ktrans]);
xlabel('x (pixels)'); ylabel('y (pixels)'); ylabel(C,'Ktrans (min-1)');

guidata(hObject, handles); % Update the GUI data structure


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.max_Ktrans = get(hObject, 'Value');

current_Ktrans_slice=handles.current_Ktrans_slice;
current_Ktrans_acquisition=handles.current_Ktrans_acquisition;

Ktrans_map = handles.Ktrans_map{current_Ktrans_acquisition,current_Ktrans_slice};
acquisition_time = handles.acquisition_time;
max_Ktrans = handles.max_Ktrans;

Ktrans_img = Ktrans_map(:,:)/acquisition_time*60; % min-1
 
axes(handles.axes2); colormap(handles.axes2,jet)
imagesc(Ktrans_img,[0 max_Ktrans]); C=colorbar; 
ylim(C,[0 max_Ktrans]);
xlabel('x (pixels)'); ylabel('y (pixels)'); ylabel(C,'Ktrans (min-1)');

guidata(hObject, handles); % Update the GUI data structure



% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DCE_repetitions_Callback(hObject, eventdata, handles)
% hObject    handle to DCE_repetitions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DCE_repetitions as text
%        str2double(get(hObject,'String')) returns contents of DCE_repetitions as a double

handles.DCE_acquisitions = str2double(get(hObject,'String'));

set(handles.current_MRI_acquisition_slider,'max',handles.DCE_acquisitions);
set(handles.slider2,'max',handles.DCE_acquisitions);

guidata(hObject, handles); % Update the GUI data structure



% --- Executes during object creation, after setting all properties.
function DCE_repetitions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DCE_repetitions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double

handles.Gd_relaxivity = str2double(get(hObject,'String'));

guidata(hObject, handles); % Update the GUI data structure


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double

handles.T1_brain_relaxation_time = str2double(get(hObject,'String'));

guidata(hObject, handles); % Update the GUI data structure


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit13_Callback(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit13 as text
%        str2double(get(hObject,'String')) returns contents of edit13 as a double

handles.T1_blood_relaxation_time = str2double(get(hObject,'String'));

guidata(hObject, handles); % Update the GUI data structure


% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double

handles.repetition_time = str2double(get(hObject,'String'));

guidata(hObject, handles); % Update the GUI data structure


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double

handles.acquisition_time = str2double(get(hObject,'String'));

guidata(hObject, handles); % Update the GUI data structure



% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function slice_number_Callback(hObject, eventdata, handles)
% hObject    handle to slice_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slice_number as text
%        str2double(get(hObject,'String')) returns contents of slice_number as a double

    handles.DCE_slices = str2double(get(hObject,'String'));
    
    set(handles.current_MRI_slice_slider,'max',handles.DCE_slices);
    set(handles.slider1,'max',handles.DCE_slices);

guidata(hObject, handles); % Update the GUI data structure



% --- Executes during object creation, after setting all properties.
function slice_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slice_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filter_matrix_size_Callback(hObject, eventdata, handles)
% hObject    handle to filter_matrix_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filter_matrix_size as text
%        str2double(get(hObject,'String')) returns contents of filter_matrix_size as a double

    handles.filter_matrix_size = str2double(get(hObject,'String'));
      
    handles.FilterMatrix = ones(handles.filter_matrix_size,handles.filter_matrix_size) / (handles.filter_matrix_size.^2);   %Smoothing pixel data with neighbors

guidata(hObject, handles); % Update the GUI data structure


% --- Executes during object creation, after setting all properties.
function filter_matrix_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filter_matrix_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

 DCE_MRI_frames = handles.DCE_MRI_frames;
 Ktrans_map = handles.Ktrans_map;
 Ktrans_map_threshold = handles.Ktrans_map_threshold;
 Kep_map = handles.Kep_map;
 t0_map = handles.t0_map;
 norm_map = handles.norm_map;
 count = handles.count;
 
 save('DCE_MRI_data.mat','DCE_MRI_frames','-mat');
 save('Ktrans_data.mat','Ktrans_map','-mat');
 save('Ktrans_thresholded_data.mat','Ktrans_map_threshold','-mat');
 save('Kep_data.mat','Kep_map','-mat');
 save('t0_data.mat','t0_map','-mat');
 save('norm_map.mat','norm_map','-magmat');
 
 guidata(hObject, handles); % Update the GUI data structure


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function current_MRI_acquisition_slider_Callback(hObject, eventdata, handles)
% hObject    handle to current_MRI_acquisition_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.current_MRI_acquisition = uint8(get(hObject, 'Value'));

current_slice = handles.current_MRI_slice;
current_acquisition = handles.current_MRI_acquisition;

axes(handles.axes1);
imagesc(handles.DCE_MRI_frames{current_acquisition,current_slice}); %colormap(gray);
%imagesc(handles.DCE_MRI_frames{current_acquisition,current_slice},[0 256]); colormap(gray);
xlabel('x (pixels)'); ylabel('y (pixels)'); C=colorbar; ylim(C,[0 256]);

delete(handles.acquisition_text);
slice_number = ['Acquisition: ',num2str(current_acquisition),'/',num2str(handles.DCE_acquisitions)];
handles.acquisition_text = annotation('textbox',[0.34 0.83 0.1 0.1],'String',slice_number,'FontSize',14,'Color','w','LineStyle','none');

guidata(hObject, handles); % Update the GUI data structure


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles.current_Ktrans_acquisition = uint8(get(hObject, 'Value'));

current_Ktrans_slice=handles.current_Ktrans_slice;
current_Ktrans_acquisition=handles.current_Ktrans_acquisition;

Ktrans_map = handles.Ktrans_map{current_Ktrans_acquisition,current_Ktrans_slice};
acquisition_time = handles.acquisition_time;
max_Ktrans = handles.max_Ktrans;

Ktrans_img = Ktrans_map(:,:)/acquisition_time*60; % min-1
 
axes(handles.axes2); colormap(handles.axes2,jet)
imagesc(Ktrans_img,[0 max_Ktrans]); C=colorbar; 
ylim(C,[0 max_Ktrans]);
xlabel('x (pixels)'); ylabel('y (pixels)'); ylabel(C,'Ktrans (min-1)');

delete(handles.acquisition_text_2);
acquisition_number_2 = ['Acquisition: ',num2str(current_Ktrans_acquisition),'/',num2str(handles.DCE_acquisitions)];
handles.acquisition_text_2 = annotation('textbox',[0.68 0.83 0.1 0.1],'String',acquisition_number_2,'FontSize',14,'Color','w','LineStyle','none');

guidata(hObject, handles); % Update the GUI data structure


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on current_MRI_acquisition_slider and none of its controls.
function current_MRI_acquisition_slider_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to current_MRI_acquisition_slider (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

set(current_MRI_acquisition_slider,'min',1);
set(current_MRI_acquisition_slider,'max',handles.DCE_acquisitions);
set(current_MRI_acquisition_slider, 'SliderStep', [1/handles.DCE_acquisitions , 5/handles.DCE_acquisitions]);

guidata(hObject, handles); % Update the GUI data structure


% --- Executes on key press with focus on current_MRI_slice_slider and none of its controls.
function current_MRI_slice_slider_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to current_MRI_slice_slider (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

set(handles.current_MRI_slice_slider,'min',1);
set(handles.current_MRI_slice_slider,'max',handles.DCE_slices);
set(handles.current_MRI_slice_slider, 'SliderStep', [1/handles.DCE_slices , 5/handles.DCE_slices]);

guidata(hObject, handles); % Update the GUI data structure


% --- Executes when selected object is changed in uibuttongroup3.
function uibuttongroup3_SelectionChangedFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uibuttongroup3 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.acq_slice_selection = get(hObject,'String');

if  strcmp(handles.acq_slice_selection,'Current slice/acquisition') == 1
    
    handles.selected_slice = handles.current_MRI_slice;
    handles.selected_acquisition = handles.current_MRI_acquisition;
    handles.selection_preference = 1;

elseif strcmp(handles.acq_slice_selection,'Selected slices/acquisitions') == 1
    
    prompt = {'Enter initial slice:','Enter final slice:','Enter reference acquisition','Enter initial acquisition','Enter final acquisition'};
    dlg_title = 'Slice / acquisition selection';
    num_lines = 1;
    defaultans = {'1','18','4','1','64'};
    selections = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
    handles.selected_initial_slice = str2num(selections{1});
    handles.selected_final_slice = str2num(selections{2});
    handles.selected_reference_acquisition = str2num(selections{3});
    handles.selected_initial_acquisition = str2num(selections{4});
    handles.selected_final_acquisition = str2num(selections{5});
   
    handles.selection_preference = 2;

end

guidata(hObject, handles); % Update the GUI data structure

    


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Ktrans_map = handles.Ktrans_map;
position = handles.ROI_position;

threshold = 0.5;

% Set a threshold to get rid of arteries %

Ktrans_threshold = zeros(180,150,handles.DCE_acquisitions);
slice = handles.current_Ktrans_slice;
    
    for acquisition = handles.selected_initial_acquisition:handles.selected_final_acquisition
                       
           Ktrans_threshold(:,:,acquisition) = Ktrans_map{acquisition,slice};
    end
    
    for i = round(position(2)):(round(position(2))+round(position(4)))
        for j = round(position(1)):(round(position(1))+round(position(3)))
            
            if Ktrans_threshold(i,j,round(handles.DCE_acquisitions/2)) > threshold*Ktrans_threshold(:,:,handles.DCE_acquisitions)
                Ktrans_threshold(i,j,:) = 0;
            end 
            
        end
    end
    

    for acquisition = handles.selected_initial_acquisition:handles.selected_final_acquisition
           handles.Ktrans_map_threshold{acquisition,slice} = Ktrans_threshold(:,:,acquisition);
    end
    
% handles.Ktrans_map = Ktrans_map;

current_Ktrans_slice=handles.current_Ktrans_slice;
current_Ktrans_acquisition=handles.current_Ktrans_acquisition;

Ktrans_map = handles.Ktrans_map_threshold{current_Ktrans_acquisition,current_Ktrans_slice};
acquisition_time = handles.acquisition_time;
max_Ktrans = handles.max_Ktrans;

Ktrans_img = Ktrans_map(:,:)/acquisition_time*60; % min-1
 
axes(handles.axes2); colormap(handles.axes2,jet);
imagesc(Ktrans_img,[0 max_Ktrans]); C=colorbar; 
ylim(C,[0 max_Ktrans]);
xlabel('x (pixels)'); ylabel('y (pixels)'); ylabel(C,'Ktrans (min-1)');

guidata(hObject, handles); % Update the GUI data structure


% --- Executes on slider movement.
function slider8_Callback(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

threshold = get(hObject,'Value');

Ktrans_map = handles.Ktrans_map;
position = handles.ROI_position;

% Set a threshold to get rid of arteries %

Ktrans_threshold = zeros(180,150,handles.DCE_acquisitions);
slice = handles.current_Ktrans_slice;
    
    for acquisition = handles.selected_initial_acquisition:handles.selected_final_acquisition
                       
           Ktrans_threshold(:,:,acquisition) = Ktrans_map{acquisition,slice};
    end
    
    for i = round(position(2)):(round(position(2))+round(position(4)))
        for j = round(position(1)):(round(position(1))+round(position(3)))
            
            if Ktrans_threshold(i,j,round(handles.DCE_acquisitions/2)) > threshold*Ktrans_threshold(:,:,handles.DCE_acquisitions)
                Ktrans_threshold(i,j,:) = 0;
            end 
            
        end
    end
    
    for acquisition = handles.selected_initial_acquisition:handles.selected_final_acquisition
                       
           handles.Ktrans_map_threshold{acquisition,slice} = Ktrans_threshold(:,:,acquisition);
           
    end
    
current_Ktrans_slice=handles.current_Ktrans_slice;
current_Ktrans_acquisition=handles.current_Ktrans_acquisition;

Ktrans_map = handles.Ktrans_map_threshold{current_Ktrans_acquisition,current_Ktrans_slice};
acquisition_time = handles.acquisition_time;
max_Ktrans = handles.max_Ktrans;

Ktrans_img = Ktrans_map(:,:)/acquisition_time*60; % min-1
 
axes(handles.axes2); colormap(handles.axes2,jet);
imagesc(Ktrans_img,[0 max_Ktrans]); C=colorbar; 
ylim(C,[0 max_Ktrans]);
xlabel('x (pixels)'); ylabel('y (pixels)'); ylabel(C,'Ktrans (min-1)');

guidata(hObject, handles); % Update the GUI data structure


% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,filepath] = uigetfile('*.mat','Select the Ktrans file');

load([filepath filename]);

current_Ktrans_slice=handles.current_Ktrans_slice;
current_Ktrans_acquisition=handles.current_Ktrans_acquisition;

Ktrans_map_to_plot = Ktrans_map{current_Ktrans_acquisition,current_Ktrans_slice};
acquisition_time = handles.acquisition_time;
max_Ktrans = handles.max_Ktrans;

Ktrans_img = Ktrans_map_to_plot(:,:)/acquisition_time*60; % min-1
 
axes(handles.axes2); colormap(handles.axes2,jet);
imagesc(Ktrans_img,[0 max_Ktrans]); C=colorbar; 
ylim(C,[0 max_Ktrans]);
xlabel('x (pixels)'); ylabel('y (pixels)'); ylabel(C,'Ktrans (min-1)');

slice_number_2 = ['Slice: ',num2str(handles.current_Ktrans_slice),'/',num2str(handles.Ktrans_slices)];
handles.slice_text_2 = annotation('textbox',[0.68 0.85 0.1 0.1],'String',slice_number_2,'FontSize',14,'Color','w','LineStyle','none');
acquisition_number_2 = ['Acquisition: ',num2str(handles.current_Ktrans_acquisition),'/',num2str(handles.Ktrans_acquisitions)];
handles.acquisition_text_2 = annotation('textbox',[0.68 0.83 0.1 0.1],'String',acquisition_number_2,'FontSize',14,'Color','w','LineStyle','none');

handles.Ktrans_map = Ktrans_map;

guidata(hObject, handles); % Update the GUI data structure
