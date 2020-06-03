clear all
cs = loadCIF('Ti-Titanium-alpha');
clc;
%% 
%The following part of the code reads abaqus .rpt files and 
%extracts orientation of the crystals as a 3x3 rotation matrix 
%and elastic strains from the simulation as a 3x3 matrix

if ~exist(mtex_path)
    error('This code requires the MTEX toolbox. (https://mtex-toolbox.github.io/')
end

folder='T:\Phani\GitHub\2Ddiffsimul\Demoset';
cd(folder)
MyFileInfo = dir('*.rpt');
[~, reindex] = sort( str2double( regexp( {MyFileInfo.name}, '\d+', 'match', 'once' )));
MyFileInfo = MyFileInfo(reindex) ;
fid = fopen(MyFileInfo(1).name);
formatSpec = '%s';
B=fscanf(fid,formatSpec);
fclose(fid);
numgrains=size(strfind(B,'Total'),2);
increments=size(MyFileInfo,1);

format long;
disp('Scanning output files');
for i=1:increments
    fid = fopen([num2str(i),'.rpt']);
    C{i} = textscan(fid,'%s','Delimiter','','EndOfLine','\r\n');
    fclose(fid);
    A=contains(C{1,i}{1,1},'Total');
    E{i}=C{1,i}{1,1}(A);
   for j=1:numgrains
       X{i,j}=strsplit(E{1,i}{j},' ');
   end
end


R=zeros(3,3,numgrains,increments);        %Rotation matrics
ee=zeros(3,3,numgrains,increments);       %Elastic strain matrics
for m=1:increments
   for k=1:numgrains
    R(1,1,k,m)=str2num(X{m,k}{1,2})/64;
    R(1,2,k,m)=str2num(X{m,k}{1,3})/64;
    R(1,3,k,m)=str2num(X{m,k}{1,4})/64;
    R(2,1,k,m)=str2num(X{m,k}{1,5})/64;
    R(2,2,k,m)=str2num(X{m,k}{1,6})/64;
    R(2,3,k,m)=str2num(X{m,k}{1,7})/64;
    R(3,1,k,m)=str2num(X{m,k}{1,8})/64;
    R(3,2,k,m)=str2num(X{m,k}{1,9})/64;
    R(3,3,k,m)=str2num(X{m,k}{1,10})/64;
    ee(1,1,k,m)=str2num(X{m,k}{1,11})/64;
    ee(1,2,k,m)=str2num(X{m,k}{1,12})/64;
    ee(1,3,k,m)=str2num(X{m,k}{1,13})/64;
    ee(2,1,k,m)=str2num(X{m,k}{1,14})/64;
    ee(2,2,k,m)=str2num(X{m,k}{1,15})/64;
    ee(2,3,k,m)=str2num(X{m,k}{1,16})/64;
    ee(3,1,k,m)=str2num(X{m,k}{1,17})/64;
    ee(3,2,k,m)=str2num(X{m,k}{1,18})/64;
    ee(3,3,k,m)=str2num(X{m,k}{1,19})/64;
   end
end
%% 


%% 
% If not using Abaqus, the above section of this code should be modified to
% input orientations of the grains as a 3x3 rotation matrix 
% and elastic strains as a 3x3 rotation matrix

% The orientation matrix R is of size 3x3xpxq, 
% where p is the number of grains and q is the number of increments
% 
% The elastic strains matrix ee of size 3x3xpxq, 
% where p is the number of grains and q is the number of increments


%hkl denote the crystal planes for which the simulation is required. 

hkl=[0,1,4;
    0,2,2;
    0,0,4;
    1,1,2;
    0,2,0;
    0,1,3;
    1,1,0;
    0,1,2;
    0,1,1;
    0,0,2;
    0,1,0];


oris=R(:,:,:,1);  %Pick the zeroth frame

rhkl=zeros([size(hkl,1),3,numgrains]);
for i=1:size(hkl,1)
    for j=1:numgrains
        rhkl(i,:,j)=oris(:,:,j)*hkl(i,:).';       
    end
end
rhkl(:,3,:)=0;
rhkl=permute(rhkl,[3,2,1]);

clear nor
for i = 1:size(rhkl,1)
    nor(i,1)=norm(rhkl(i,:,1));
end
nor=repmat(nor,[1,3,size(hkl,1)]);
%s0=zeros([size(hkl,1),3,noofori]);

s0= (1./nor).*rhkl;
%omega=atan(s0(:,2,:)./s0(:,1,:));
%omega1=pi/2+(pi/2-omega);
%omega=cat(1,omega,pi/2+(pi/2-omega));

%omega=repmat(omega,[1 increments 1]);
s0=repmat(s0,[2,1,1]);
s0(numgrains+1:2*numgrains,1,:)=s0(numgrains+1:2*numgrains,1,:).*-1;
omega=atan(s0(:,2,:)./s0(:,1,:));
omega(numgrains+1:2*numgrains,1,:)=pi/2+(pi/2-omega(numgrains+1:2*numgrains,1,:));

d0spacing=zeros([2*numgrains,increments,size(hkl,1)]);
    for i =1:size(hkl,1)
        m1=Miller(hkl(i,1),hkl(i,2),hkl(i,3),cs);
        d0spacing(:,1,i)=m1.dspacing;
    end
    
dspacing=zeros(2*numgrains,increments,size(hkl,1));    
for j = 1:2*numgrains
    for i= 1:increments
        for k = 1:size(hkl,1)
        ee_frame=ee(:,:,:,i);
        ee_frame=repmat(ee_frame,[1,1,2]);
        dspacing(j,i,k)=(1+s0(j,:,k)*ee_frame(:,:,j)*s0(j,:,k)')*d0spacing(j,1,k);
        end
    end
end



fig = figure;
for i = 1:size(hkl,1)
    polarplot(omega(:,1,i),d0spacing(:,1,i),'.','Color','w','MarkerSize',4);
    %[x,y]=pol2cart(omega(:,1,i),repmat(dspacing,[1,numgrains])');
    %scatter(x,y);
    hold on
    %polarplot(omega1(:,1,i),repmat(dspacing,[1,numgrains]),'.','Color','w','MarkerSize',4);
end
ax = gca;
ax.RGrid = 'off';
ax.ThetaGrid = 'off';
ax.ThetaAxis.Visible='off';
ax.RAxis.Visible='off';
ax.Color = 'k';
hold off
%fig = gcf;
fig.PaperUnits = 'points';
fig.PaperPosition = [0 0 2048 2048];
fig.Color='k';
print(['image0'],'-dpng','-r0')
close all


for j = 1:increments
    fig = figure;
    for i = 1:size(hkl,1)
        polarplot(omega(:,1,i),dspacing(:,j,i),'.','Color','w','MarkerSize',4);
        %[x,y]=pol2cart(omega(:,1,i),repmat(dspacing,[1,numgrains])');
        %scatter(x,y);
        hold on
        %polarplot(omega1(:,1,i),repmat(dspacing,[1,numgrains]),'.','Color','w','MarkerSize',4);
    end
    ax = gca;
    ax.RGrid = 'off';
    ax.ThetaGrid = 'off';
    ax.ThetaAxis.Visible='off';
    ax.RAxis.Visible='off';
    ax.Color = 'k';
    hold off
    %fig = gcf;
    fig.PaperUnits = 'points';
    fig.PaperPosition = [0 0 2048 2048];
    fig.Color='k';
    print(['image' int2str(j)],'-dpng','-r0')
    close all
end
%% 





