% this script read a *.p_2d file and plot the velocity field

clear all;
close all;

%% Initialize variables.
[FileName,PathName,FilterIndex] = uigetfile('*.p_2d');

filename = strcat(PathName,FileName);
delimiter = ' ';

%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, ...
    'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
Sim170001 = [dataArray{1:end-1}];
%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

x = Sim170001(:,1);
temp = find(x(2:end)==x(1));
nx = temp(1);
x = reshape(x(1:end-2),nx,[]);
x = x(1:end-1,:);

y = Sim170001(:,2);
y = reshape(y(1:end-2),nx,[]);
y = y(1:end-1,:);

h = Sim170001(:,3);
h = reshape(h(1:end-2),nx,[]);
h = h(1:end-1,:);

u = Sim170001(:,4);
u = reshape(u(1:end-2),nx,[]);
u = u(1:end-1,:);

v = Sim170001(:,5);
v = reshape(v(1:end-2),nx,[]);
v = v(1:end-1,:);

B = Sim170001(:,6);
B = reshape(B(1:end-2),nx,[]);
B = B(1:end-1,:);

hB = Sim170001(:,7);
hB = reshape(hB(1:end-2),nx,[]);
hB = hB(1:end-1,:);

%%
figure

mod_vel = sqrt( u.^2 + v.^2 );

surf(x,y,B,mod_vel,'EdgeColor','none');

shading interp;

colorbar;
  
box on;
 
axis equal;
view(52,36);

light;
light;

dx = x(2,1)-x(1,1);
dy = y(1,2)-y(1,1);

%%
figure;

quiver(x,y,u,v);
axis equal;


