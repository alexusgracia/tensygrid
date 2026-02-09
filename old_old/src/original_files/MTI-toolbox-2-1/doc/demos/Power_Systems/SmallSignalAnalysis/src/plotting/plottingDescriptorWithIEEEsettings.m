% 
% Plotting the descriptor figures of the paper according to IEEE standard
% KauChr
% 01.05.2025


%
%%
% According to the standard of IEEE Transactions and Journals: 

% Times New Roman is the suggested font in labels. 

% For a singlepart figure, labels should be in 8 to 10 points,
% whereas for a multipart figure, labels should be in 8 points.

% Width: column width: 8.8 cm; page width: 18.1 cm.

%%  

descriptorVload=sqrt((descriptorSimout.y(:,indexVDload+MTImodel.n).^2).'+(descriptorSimout.y(:,indexVDload+MTImodel.n+1).^2).')./Vpeak;

%GFM
descriptorGFMcurrentAmplitude=sqrt((descriptorSimout.y(:,indexGFMiD).^2).'+(descriptorSimout.y(:,indexGFMiD+1).^2).')./Ibase;
descriptorGFMpq=descriptorSimout.y(:,MTImodel.n+(indexGFMp:indexGFMp+1))./Sbase;
descriptorGFMfrequency=descriptorSimout.y(:,MTImodel.n+indexGFMomega)/(2*pi);
descriptorGFMangle=(descriptorSimout.y(:,indexGFMangle)-descriptorSimout.y(:,indexGRIDangle))*180/pi;

descriptorGFLcurrentAmplitude=sqrt((descriptorSimout.y(:,indexGFLiD).^2).'+(descriptorSimout.y(:,indexGFLiD+1).^2).')./Ibase;
descriptorGFLpq=descriptorSimout.y(:,(indexGFLp:indexGFLp+1)+MTImodel.n)./Sbase;
descriptorGFLfrequency=descriptorSimout.y(:,MTImodel.n+indexGFLomega)/(2*pi);
descriptorGFLangle=(descriptorSimout.y(:,indexGFLangle)-descriptorSimout.y(:,indexGRIDangle))*180/pi;

%% Plotting settings

% window size
 textwidth=516/2;
 conv_pt2cm=0.0351;
 wd=textwidth*conv_pt2cm;
 wdIEEE=8.8;
 widthIEEEpt=wdIEEE/conv_pt2cm;
 heightIEEEpt=696*0.4;
 h=696*conv_pt2cm/3;

% positioning on the screen
 pw=3*wd;   
 ph=5;   

% fontsize
fs=8;
tUP=15; % upper limit of the time for the xaxis in second


%% figure margins
top = 0.5;  % normalized top margin
bottom = 3;	% normalized bottom margin
left = 3.5;	% normalized left margin
right = 1;  % normalized right margin

%%



%% 
lowerTimeLimit=descriptorSimout.tsim(1);
upperTimeLimit=descriptorSimout.tsim(end);


%% --- GFM Figure ---
% Define subplot positions for 6 plots, last one higher for legend
% Reserve space for sgtitle by shifting subplots down
nSubplots = 6;
subplotGap = 0.035;
subplotGapPos = 0.045;
bottomMargin = 0.0;
lastSubplotHeight = 0.2; % Make last subplot higher for legend
subplotHeight = 0.12;

sgtitleSpace = 0.08; % normalized units, space reserved at top for sgtitle
subplotHeight = (1 - sgtitleSpace - subplotGap * (nSubplots) - bottomMargin - lastSubplotHeight) / (nSubplots - 1);

width = 0.8;
leftMargin = 0.15;

subplotPositions = cell(1, nSubplots);
for i = 1:(nSubplots-1)
    % positioning from top to bottom, shifted down by sgtitleSpace
    subplotPositions{i} = [leftMargin, ...
        1 - sgtitleSpace - (i) * subplotHeight - (i-1) * subplotGapPos + 0.01, ...
        width, subplotHeight];
end
subplotPositions{nSubplots} = [leftMargin, bottomMargin, width, lastSubplotHeight];
%% GFM figure

figure%(4); clf;
f = gcf;
% Increase figure height to compensate for sgtitle space
sgtitleSpace = 0.08; % must match the value above
heightIEEEptWithTitle = heightIEEEpt / (1 - sgtitleSpace);
set(f,'Units', 'points', 'Position', [100, 100, widthIEEEpt, heightIEEEptWithTitle],'Name','GFM NTI MTI LTI');
set(f,'defaultLineLineWidth',1.2);
set(f,'defaultAxesFontName','Times New Roman');
set(f,'defaultAxesFontSize',8);

% Add a super title above all subplots, now with reserved space
sgtitle('Grid-forming Converter','FontSize',10,'FontWeight','bold')
set(f,'defaultTextFontName','Times New Roman');
set(f,'defaultTextFontSize',8);
set(f,'defaultLegendFontName','Times New Roman');
set(f,'defaultLegendFontSize',8);

% --- P only subplot ---
ax1 = subplot('Position', subplotPositions{1});
grid on; hold on; box on;
plot(timeSimulink, gfmPQSimulink(:,1))
plot(timeDMSIM, gfmPQDMSIM(:,1), '--')
plot(descriptorSimout.tsim, descriptorGFMpq(:,1), ':','LineWidth',1.6)
ylabel('$p$ (p.u.)', 'Interpreter', 'latex')
xlim([lowerTimeLimit upperTimeLimit])
xticklabels(ax1, {});

% --- Q only subplot ---
axQ = subplot('Position', subplotPositions{2});
grid on; hold on; box on;
plot(timeSimulink, gfmPQSimulink(:,2))
plot(timeDMSIM, gfmPQDMSIM(:,2), '--')
plot(descriptorSimout.tsim, descriptorGFMpq(:,2), ':','LineWidth',1.6)
ylabel('$q$ (p.u.)', 'Interpreter', 'latex')
xlim([lowerTimeLimit upperTimeLimit])
xticklabels(axQ, {});

% --- Voltage subplot ---
ax2 = subplot('Position', subplotPositions{3});
grid on; hold on; box on;
plot(timeSimulink, voltagePCCAmplitudeSimulink)
plot(timeDMSIM, voltageAmplitudeDMSIM, '--')
plot(descriptorSimout.tsim, descriptorVload, ':','LineWidth',1.6)
ylabel('$||\mathbf{v}_{DQ}||$ (p.u.)', 'Interpreter', 'latex')
xlim([lowerTimeLimit upperTimeLimit])
xticklabels(ax2, {});

% --- Current subplot ---
ax3 = subplot('Position', subplotPositions{4});
grid on; hold on; box on;
plot(timeSimulink, gfmCurrentAmplitudeSimulink)
plot(timeDMSIM, gfmCurrentAmplitudeDMSIM, '--')
plot(descriptorSimout.tsim, descriptorGFMcurrentAmplitude, ':','LineWidth',1.6)
ylabel('$||\mathbf{i}_{DQ}||$ (p.u.)', 'Interpreter', 'latex')
xlim([lowerTimeLimit upperTimeLimit])
xticklabels(ax3, {});

% --- Frequency subplot ---
ax4 = subplot('Position', subplotPositions{5});
grid on; hold on; box on;
plot(timeSimulink, frequencySimulink)
plot(timeDMSIM, gfmFrequencyDMSIM, '--')
plot(descriptorSimout.tsim, descriptorGFMfrequency, ':','LineWidth',1.6)
ylabel('$f$ (Hz)', 'Interpreter', 'latex')
xlim([lowerTimeLimit upperTimeLimit])
xticklabels(ax4, {});

% --- Angle subplot (higher for legend) ---
ax5 = subplot('Position', subplotPositions{6});
grid on; hold on; box on;
plot(timeSimulink, deltaAngleSimulink)
plot(timeDMSIM, deltaAngleDMSIM, '--')
plot(descriptorSimout.tsim, descriptorGFMangle, ':','LineWidth',1.6)
ylabel('$\delta (\circ)$', 'Interpreter', 'latex')
xlim([lowerTimeLimit upperTimeLimit])

xlabel(ax5, '$t$ (s)', 'FontSize', 8, 'Interpreter', 'latex')
legend('Nonlinear','Multilinear','Linear','Location','southoutside','Orientation','horizontal','NumColumns',3)




%% --- GFL Figure ---

figure%(5)%; clf;
f = gcf;
% Increase figure height to compensate for sgtitle space (must match GFM logic)
sgtitleSpace = 0.08; % must match the value above
heightIEEEptWithTitle = heightIEEEpt / (1 - sgtitleSpace);
set(f,'Units', 'points', 'Position', [100, 100, widthIEEEpt, heightIEEEptWithTitle],'Name','GFL NTI MTI LTI');
set(f,'defaultLineLineWidth',1.2);
set(f,'defaultAxesFontName','Times New Roman');
set(f,'defaultAxesFontSize',8);
set(f,'defaultTextFontName','Times New Roman');
set(f,'defaultTextFontSize',8);
set(f,'defaultLegendFontName','Times New Roman');
set(f,'defaultLegendFontSize',8);
% Add a super title above all subplots, now with reserved space
sgtitle('Grid-following Converter','FontSize',10,'FontWeight','bold')
% Use the same subplotPositions as above for equal spacing
% --- P only subplot ---
ax1 = subplot('Position', subplotPositions{1});
grid on; hold on; box on;
plot(timeSimulink, gflPQSimulink(:,1))
plot(timeDMSIM, gflPQDMSIM(:,1), '--')
plot(descriptorSimout.tsim, descriptorGFLpq(:,1), ':','LineWidth',1.6)
ylabel('$p$ (p.u.)', 'Interpreter', 'latex')
xlim([lowerTimeLimit upperTimeLimit])
xticklabels(ax1, {});

% --- Q only subplot ---
axQ = subplot('Position', subplotPositions{2});
grid on; hold on; box on;
plot(timeSimulink, gflPQSimulink(:,2))
plot(timeDMSIM, gflPQDMSIM(:,2), '--')
plot(descriptorSimout.tsim, descriptorGFLpq(:,2), ':','LineWidth',1.6)
ylabel('$q$ (p.u.)', 'Interpreter', 'latex')
xlim([lowerTimeLimit upperTimeLimit])
xticklabels(axQ, {});

% --- Voltage subplot ---
ax2 = subplot('Position', subplotPositions{3});
grid on; hold on; box on;
plot(timeSimulink, voltagePCCAmplitudeSimulink)
plot(timeDMSIM, voltageAmplitudeDMSIM, '--')
plot(descriptorSimout.tsim, descriptorVload, ':','LineWidth',1.6)
ylabel('$||\mathbf{v}_{DQ}||$ (p.u.)', 'Interpreter', 'latex')
xlim([lowerTimeLimit upperTimeLimit])
xticklabels(ax2, {});

% --- Current subplot ---
ax3 = subplot('Position', subplotPositions{4});
grid on; hold on; box on;
plot(timeSimulink, gflCurrentAmplitudeSimulink)
plot(timeDMSIM, gflCurrentAmplitudeDMSIM, '--')
plot(descriptorSimout.tsim, descriptorGFLcurrentAmplitude, ':','LineWidth',1.6)
ylabel('$||\mathbf{i}_{DQ}||$ (p.u.)', 'Interpreter', 'latex')
xlim([lowerTimeLimit upperTimeLimit])
xticklabels(ax3, {});

% --- Frequency subplot ---
ax4 = subplot('Position', subplotPositions{5});
grid on; hold on; box on;
plot(timeSimulink, gflFrequencySimulink)
plot(timeDMSIM, gflFrequencyDMSIM, '--')
plot(descriptorSimout.tsim, descriptorGFLfrequency, ':','LineWidth',1.6)
ylabel('$f$ (Hz)', 'Interpreter', 'latex')
xlim([lowerTimeLimit upperTimeLimit])
xticklabels(ax4, {});

% --- Angle subplot (higher for legend) ---
ax5 = subplot('Position', subplotPositions{6});
grid on; hold on; box on;
plot(timeSimulink, gfldeltaAngleSimulink)
plot(timeDMSIM, gfldeltaAngleDMSIM, '--')
plot(descriptorSimout.tsim, descriptorGFLangle, ':','LineWidth',1.6)
ylabel('$\delta (\circ)$', 'Interpreter', 'latex')
xlim([lowerTimeLimit upperTimeLimit])
xlabel(ax5, '$t$ (s)', 'FontSize', 8, 'Interpreter', 'latex')
legend('Nonlinear','Multilinear','Linear','Location','southoutside','Orientation','horizontal','NumColumns',3)




%% 


%% Bifurcation plot

if contains(SimulationSetup.caseStudy,'Bifurcation')


    figure%(30); clf;
    f = gcf;
    set(f,'Units', 'points', 'Position', [100, 100, widthIEEEpt, 0.3*heightIEEEpt],'Name','bifurcation');
    set(f,'defaultLineLineWidth',1.2);
    set(f,'defaultAxesFontName','Times New Roman');
    set(f,'defaultAxesFontSize',8);
    set(f,'defaultTextFontName','Times New Roman');
    set(f,'defaultTextFontSize',8);
    set(f,'defaultLegendFontName','Times New Roman');
    set(f,'defaultLegendFontSize',8);
    
    % --- P only subplot ---
    %ax1 = subplot('Position', subplotPositions{1});
    grid on; hold on; box on;
    plot(timeSimulink, gfmPQSimulink(:,1))
    plot(timeDMSIM, gfmPQDMSIM(:,1), '--')
    % plot(descriptorSimout.tsim, descriptorGFMpq(:,1), ':')
     ylabel('$p$ (p.u.)', 'Interpreter', 'latex')
     xlim([lowerTimeLimit upperTimeLimit])
    % xticklabels(ax1, {});
    
    % --- Q only subplot ---
    %axQ = subplot('Position', subplotPositions{2});
    % grid on; hold on; box on;
    % plot(timeSimulink, gfmPQSimulink(:,2))
    % plot(timeDMSIM, gfmPQDMSIM(:,2), '--')
    % % plot(descriptorSimout.tsim, descriptorGFMpq(:,2), ':')
    % ylabel('$Q$ (p.u.)', 'Interpreter', 'latex')
    % xlim([lowerTimeLimit upperTimeLimit])
    % xticklabels(axQ, {});
    
    
    xlabel( '$t$ (s)', 'FontSize', 8, 'Interpreter', 'latex')
    legend('Nonlinear','Multilinear','Location','southoutside','Orientation','horizontal','NumColumns',2)


end