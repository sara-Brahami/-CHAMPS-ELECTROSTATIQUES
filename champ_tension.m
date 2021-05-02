clear variables;
 close all
clc;

fName = 'tension.txt';

%% Read data 1 potentiel
%-- Format string for each line of text
formatSpec = '%f';            delimiter = ',';
%-- Open the text file
fileID = fopen( fName, 'r' );
%-- Read columns of data according to format string.
dataArray = textscan( fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue', NaN, 'ReturnOnError', false );
%-- Close the text file
fclose( fileID );
%-- Allocate imported array to column variable names
x = dataArray{:, 1}; 
%-- Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;
%% read data 2 champ
fName = 'champs.txt';
%-- Format string for each line of text
formatSpec = '%f';            delimiter = ',';
%-- Open the text file
fileID = fopen( fName, 'r' );
%-- Read columns of data according to format string.
dataArray = textscan( fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue', NaN, 'ReturnOnError', false );
%-- Close the text file
fclose( fileID );
%-- Allocate imported array to column variable names
z = dataArray{:, 1}; 
%-- Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

%% converting
sizee=sqrt(length(x));
v=zeros(sizee,sizee);
for i=0:sizee-1
v(i+1,:)=x([sizee*i+1:sizee*(i+1)]);
end
%% converting 2
sizee1=sqrt(length(z)/2);
Ex=zeros(sizee,sizee);
Ey=zeros(sizee,sizee);
for i=0:sizee-1
Ex(i+1,:)=z([sizee*i+1:sizee*(i+1)]);
Ey(i+1,:)=z([sizee^2+sizee*i+1:sizee^2+sizee*(i+1)]);
end

%% plot
figure(1);
contour(v,20);colorbar;title("Le potentiel et le champs en 2D ");
hold on;quiver(Ex,Ey);
figure(2);
surf(v);colorbar;title("Le potentiel en 3D ");