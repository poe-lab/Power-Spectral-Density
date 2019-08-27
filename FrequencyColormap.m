function varargout = FrequencyColormap(varargin)
% FREQUENCYCOLORMAP M-file for FrequencyColormap.fig
%      FREQUENCYCOLORMAP, by itself, creates a new FREQUENCYCOLORMAP or raises the existing
%      singleton*.
%
%      H = FREQUENCYCOLORMAP returns the handle to a new FREQUENCYCOLORMAP or the handle to
%      the existing singleton*.
%
%      FREQUENCYCOLORMAP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FREQUENCYCOLORMAP.M with the given input arguments.
%
%      FREQUENCYCOLORMAP('Property','Value',...) creates a new FREQUENCYCOLORMAP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FrequencyColormap_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FrequencyColormap_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FrequencyColormap

% Last Modified by GUIDE v2.5 27-Sep-2010 15:52:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FrequencyColormap_OpeningFcn, ...
                   'gui_OutputFcn',  @FrequencyColormap_OutputFcn, ...
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


% --- Executes just before FrequencyColormap is made visible.
function FrequencyColormap_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FrequencyColormap (see VARARGIN)

% Choose default command line output for FrequencyColormap
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FrequencyColormap wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FrequencyColormap_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in browseButton.
function browseButton_Callback(hObject, eventdata, handles)
global matFile
working_dir=pwd;
current_dir='C:\SleepData\Results';
cd(current_dir);
[filename, pathname] = uigetfile('*.mat', 'Select the place field file.');
if isequal(filename,0) || isequal(pathname,0)
    uiwait(errordlg('You need to select a file. Please try again',...
        'ERROR','modal'));
    cd(working_dir);
else
    cd(working_dir);
    matFile= fullfile(pathname, filename);
    set(handles.matFileName,'string',filename);
end


% --- Executes on button press in analyzeButton.
function analyzeButton_Callback(hObject, eventdata, handles)
global matFile
load(matFile, 'Occ', 'TC');
CellOccupancyMatrix = Occ; % A matrix of occupancy calculations for each bin of the maze
cellTotalSpikeMatrix = TC;
clear Occ TC
numberOfCells = size(cellTotalSpikeMatrix, 2); % Code for finding length of array
gridSize = size(CellOccupancyMatrix);
numberRows = gridSize(1);
numberColumns = gridSize(2);

for i = 1:numberOfCells    
    normFreq = cellTotalSpikeMatrix{i}./CellOccupancyMatrix;
    for j=1:numberRows
        for k=1:numberColumns
            if isequalwithequalnans(normFreq(j,k),NaN)
                if isequal(CellOccupancyMatrix(j,k), 0)
                    CellOccupancyMatrixB(j,k) = NaN;
                    normFreq(j,k) = NaN;
                
                end 
            end
        end
    end
    figure(i)
    axes('FontSize',20,'FontName','Arial')
    surf(normFreq)
    axis([0 numberRows 0 numberColumns 0 150 0 20])
    colorbar('location','eastoutside')
    xlabel('Pixels')
    ylabel('Pixels')
    zlabel('Hz')
    title(['Cell ' num2str(i) ' Average Firing Rate'])
end
 for j=1:numberRows
    for k=1:numberColumns
        if isequal(CellOccupancyMatrix(j,k), 0)
            CellOccupancyMatrix(j,k) = NaN;
        end
    end
 end
                                        
figure(i+1)
axes('FontSize',20,'FontName','Arial')
surf(CellOccupancyMatrix)
axis([0 numberRows 0 numberColumns])
colorbar('location','eastoutside')
xlabel('Pixels')
ylabel('Pixels')
zlabel('Hz')
title(['Occupancy Colormap (s)'])
