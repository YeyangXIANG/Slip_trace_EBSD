%% Slip trace analysis

%%

clear,clc;close all

%% Define parameters

% Number of strace to be analyzed
num_strace = 1;

% Phase to be analyzed
Phase_chosen = 'Magnesium';

% Define x and y directions to match the SEM images
plotx2east; plotzIntoPlane

% Degree threshold of grain boundary
GrainAngle = 10;

%% Load EBSD data

[fname, pname] = uigetfile('*.crc', 'Choose ebsd .crc');
cd(pname)

ebsd = EBSD.load([pname '\' fname],'convertEuler2SpatialReferenceFrame','setting 2');
ebsd_raw = ebsd;
CS = ebsd.CS;

% plot raw map
figure
[~,mP] = plot(ebsd,ebsd.orientations);
print(gcf,[fname(1:end-4) '_RawMap'],'-dpng','-r400');
% saveas(gcf,[fname(1:end-4) '_RawMap'],'fig');

mP.micronBar.visible = 'off';
print(gcf,[fname(1:end-4) '_RawMap_NoScaleBar'],'-dpng','-r400');

% plot ipf key
ipfKey = ipfColorKey(ebsd(Phase_chosen));
figure;
plot(ipfKey)
print(gcf,[fname(1:end-4) '_ipfKey'],'-dpng','-r400');

%% Reconstruct grains

%%%% Plot raw grains map
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',GrainAngle*degree);
ebsd(grains(grains.grainSize<=10)) = [];
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd('indexed'),'angle',GrainAngle*degree);

figure;
plot(ebsd,ebsd.orientations);
hold on
plot(grains.boundary,'linewidth',1);
hold off
saveas(gcf,[fname(1:end-4) '_RawGrains'],'fig');
print(gcf,[fname(1:end-4) '_RawGrains'],'-dpng','-r400');

%%%% fill missing
F = halfQuadraticFilter; %F = infimalConvolutionFilter; %F = splineFilter;
ebsd = smooth(ebsd('indexed'),F,'fill');     % fill all
% ebsd = smooth(ebsd,F,'fill',grains);              % only use indexded part

% redo grains reconstruction
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',GrainAngle*degree);
% smooth
grains = grains.smooth(5);

% extract new grain boundaries
gB = grains.boundary;
gB_MgMg = gB(Phase_chosen,Phase_chosen);

% plot denoised grains

figure;
plot(ebsd,ebsd.orientations);
hold on
plot(grains.boundary,'linewidth',1);
hold off
% saveas(gcf,[fname(1:end-4) '_DenoiseMap'],'fig');
print(gcf,[fname(1:end-4) '_DenoiseMap'],'-dpng','-r400');

%% Load corresponding SEM image 
% from this section, run many times for different SEM images and save them together

warning off
[ImagName, ImagPath] = uigetfile('*.tif', 'Choose SEM image');
cd(ImagPath);
n = 3;

%% Slip Trace Analysis

start = 1;
over = num_strace;
% For each slip trace, store 3 most likely slip system
for i = start:over
    nn = 1;
    while nn == 1
        [activePlanes,traceMismatch,grainID,line] = SlipTraceHCP(ImagName,ebsd,grains,n)
        % The choise of optional analysis 
        prompt = ['Try again for this grain? 1:Yes, 2:No'];
        titlename = ['Fit slip trace'];
        nn = inputdlg(prompt,titlename);
        nn = str2num(nn{1});
        close(100);
        close(101);
        close(102);
        close(103);
    end
    PlaneAll{i}(1,:) = activePlanes.h;
    PlaneAll{i}(2,:) = activePlanes.k;
    PlaneAll{i}(3,:) = activePlanes.i;
    PlaneAll{i}(4,:) = activePlanes.l;
    grainIDALL(i) = grainID;
    traceMismatchAll{i} = traceMismatch;
    lineAll{i} = line;
end

%%

% cs = crystalShape.hex(ebsd.CS);
% 
% color = ipfKey.orientation2color(ebsd(grains(grainIDALL(i))).orientations);
% ipfKey = ipfColorKey(ebsd(Phase_chosen));
% % set the referece direction to X
% ipfKey.inversePoleFigureDirection = vector3d.X;
% % compute the colors
% color = ipfKey.orientation2color(ebsd(grains(grainIDALL(i))).orientations);
% 
% figure
% plot(ebsd(grains(grainIDALL(i))),color,'micronbar','off');
% hold on
% plot(grains(grainIDALL(i)).boundary,'linewidth',2)
% plot(grains(grainIDALL(i)),0.5*cs,'faceAlpha',0.3,'linewidth',2);
% hold off
% box off
% axis off
% print(gcf,[ImagName(1:end-4) '_grain' num2str(grainID) 'orientation'],'-dpng','-r400');

%% Extract the best fit

% Extract data
for i = 1:length(PlaneAll)
    PlaneBest{i} = PlaneAll{i}(:,1);
    % Replace specific values by 0
    PlaneBest{i}(find(abs(PlaneBest{i}) < 0.001)) = 0;
end

% Convert to string
for i = 1:length(PlaneAll)
    PlaneString{i} = '[';
    for j = 1:4
        PlaneString{i} = append(PlaneString{i},num2str(PlaneBest{i}(j,1)));
    end
    PlaneString{i} = append(PlaneString{i},']');
end

%% Plot results on ebsd map
figure;
% plot(ebsd,ebsd.orientations)
% hold on
% plot(grains.boundary,'linewidth',1)
% hold off
plot(grains,grains.meanOrientation);
cs = crystalShape.hex(ebsd.CS);
for i = 1:length(PlaneAll)
    hold on
    plot(grains(grainIDALL(i)),0.5*cs,'faceAlpha',0.1);
    text(median(grains(grainIDALL(i)).x),median(grains(grainIDALL(i)).y),...
        ['Grain',num2str(i),' ', PlaneString{i}],'fontsize',13,'FontWeight','bold')
end
hold off
print(gcf,[ImagName(1:end-4) '_' num2str(i)],'-dpng','-r400');
warning on

%% Save results

savename = [ImagName(1:end-4) '_SlipTrace_Analyzed'];
save(savename,'ebsd','grains','grainIDALL','PlaneAll','PlaneBest','PlaneString','traceMismatchAll')

%% pole figure for specific grain
% ebsd_grain = ebsd(grains(150));
% % choose the plane for mapping
% h = [Miller(0,0,0,1,ebsd_grain(Phase_chosen).CS),...
%     Miller(1,0,-1,0,ebsd_grain(Phase_chosen).CS),...
%     Miller(1,0,-1,1,ebsd_grain(Phase_chosen).CS),...
%     Miller(1,1,-2,2,ebsd_grain(Phase_chosen).CS)];%,Miller(1,0,-1,2,ebsd(Phase_chosen).CS)];
% 
% % plot pole figure
% figure;
% plotPDF(ebsd_grain(Phase_chosen).orientations,h,'contourf','minmax','antipodal');
% % add colorbar
% CLim(gcm,'equal');
% mtexColorbar

%% Local function

function [activePlanes,traceMismatch,grainID,line] = SlipTraceHCP(ImagName,ebsd,grains,n)

% Version 2.0, 2021/04/23, Yeyang
% get the n most possible slip planes and the corresponding mismatch angles
% return the chosen slip trace line, and the ID of grain chosen

% % To run the function:
% read image file and zoom in
% click two end points of slip trace
% choose corresponding grain

%%

if(nargin<4)
    n=3;
end

%% Extract slip trace
%Read in image file and have user click on slip trace
image=imread(ImagName);
figure(100);
imshow(image);

% choose a part to zoom in
%uiwait(msgbox('Click the lower left point and then upper right point to zoom in','Select Corners','modal'));
title('Click the lower left point and then upper right point to zoom in');
P = get(gca);
P.OuterPosition(4)=P.OuterPosition(4)*0.98;
set(gca,'OuterPosition',P.OuterPosition)
scrsz = get(0,'ScreenSize');
set(gcf,'Position',scrsz);
zoom1 = ginput(1);
zoom2 = ginput(1);
figure(101);
imshow(image(zoom2(2):zoom1(2),zoom1(1):zoom2(1)))

%uiwait(msgbox('Click two points on a single slip trace.','Select Corners','modal'));
title('Click two points on a single slip trace')
P2 = get(gca);
P2.OuterPosition(4)=P2.OuterPosition(4)*0.98;
set(gca,'OuterPosition',P2.OuterPosition)
scrsz = get(0,'ScreenSize');
set(gcf,'Position',scrsz);
c1 = ginput(1);
c2 = ginput(1);
line=normalize(vector3d(c2(1)-c1(1),c2(2)-c1(2),0));

%% Choose corresponding grain

CS = grains.CS;

%Display plot of grains and have user select the one that their SEM image is located in.
plotx2east; plotzIntoPlane;%( x to right and y to down, corresponds to SEM line)

%%%% Option 1, use grain mean orientation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(102)
% plot(grains,grains.meanOrientation);
% g=ginput(1);
% 
% grainID = grains.findByLocation(g);
% grain=grains(grainID);
% ori=grain.meanOrientation;

%%%% Option 2, use the local orientation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(102)
plot(ebsd,ebsd.orientations);
hold on
plot(grains.boundary,'linewidth',1.5);
hold off
g=ginput(1);
grainID = grains.findByLocation(g);
ori = ebsd(ebsd.findByLocation(g)).orientations;

%% Define slip systems

Basal = slipSystem.basal(CS).symmetrise('antipodal');
Basal = Basal(1);
Basal_r = ori * Basal;

PrismaticA = slipSystem.prismaticA(CS).symmetrise('antipodal');
PrismaticA_r = ori * PrismaticA;

% % It seems we can't tell the difference between PyramidalA and PyramidalCA by slip trace
% PyramidalA = slipSystem.pyramidalA(CS).symmetrise('antipodal');
% PyramidalA = ori * PyramidalA;

PyramidalCA = slipSystem.pyramidalCA(CS).symmetrise('antipodal');
PyramidalCA = PyramidalCA([1,3,5,7,9,11]);
PyramidalCA_r = ori * PyramidalCA;

Pyramidal2CA = slipSystem.pyramidal2CA(CS).symmetrise('antipodal');
Pyramidal2CA_r = ori * Pyramidal2CA;

%% Plot trace of slip systems

figure(103)
for i = 1:length(Basal_r)
    hold on
    h(1) = quiver(grains(grainID),Basal_r(i).trace,'color','g');
end
for i = 1:length(PrismaticA_r)
    hold on
    h(2) = quiver(grains(grainID),PrismaticA_r(i).trace,'color','b');
end
for i = 1:length(PyramidalCA_r)
    hold on
    h(3) = quiver(grains(grainID),PyramidalCA_r(i).trace,'color','c');
end
for i = 1:length(Pyramidal2CA_r)
    hold on
    h(4) = quiver(grains(grainID),Pyramidal2CA_r(i).trace,'color',[0.9290 0.6940 0.1250]);
end
hold on
h(5) = quiver(grains(grainID),line,'color','k');
quiver(grains(grainID),-line,'color','k');
hold off

set(gca,'YDir','reverse');  
axis equal
axis off

savename = [ImagName(1:end-4) '_grain' num2str(grainID) '_SlipTrace'];
print(gcf,savename,'-dpng','-r400');

% save image
legend(h,{'Basal';'Prismatic';'Pyramidal I';'Pyramidal II';'Slip trace'},'edgecolor','none','location','northeast')

%% Get the minimum angle and the slip plane

ang_basal = angle(line,Basal_r.trace);
ang_prismatic = angle(line,PrismaticA_r.trace);
ang_pyramidalCA = angle(line,PyramidalCA_r.trace);
ang_pyramidal2CA = angle(line,Pyramidal2CA_r.trace);

% All angles, 1;basal, 2-4:prismatic, 5-10:pyramidalI, 11-16:pyramidalII
angles = [ang_basal,ang_prismatic,ang_pyramidalCA,ang_pyramidal2CA];

planeAll = [Basal.n; PrismaticA.n; PyramidalCA.n; Pyramidal2CA.n];

% Get the n slip planes that produce traces closest to the trace clicked on
[ASorted AIdx] = sort(angles);
if n<length(AIdx)
    top = AIdx(1:n);
else
    top=AIdx;
end

activePlanes = planeAll(top);
traceMismatch = angles(top)/degree;

text(0.1,0.1,['Mismatch degree = ' num2str(traceMismatch)],'Units','normalized','Color','red','FontSize',12);

end
