% 
% Plotting the figures of the paper according to IEEE standard
% KauChr
% 01.05.2025
%


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

set(0,'defaultLineLineWidth',1.2);

%% grid-forming converter

                
              
                figure%(1); clf;
                f = gcf;
                % Reserve space for sgtitle and increase figure height
                sgtitleSpace = 0.08;
                heightIEEEptWithTitle = heightIEEEpt / (1 - sgtitleSpace);
                set(f,'Units', 'points', 'Position', [100, 100, widthIEEEpt, heightIEEEptWithTitle],'Name','GFM Simulink MTI');

                set(f,'defaultLineLineWidth',1.2);
                
                set(f,'defaultAxesFontName','Times New Roman');
                set(f,'defaultAxesFontSize',8);
                
                set(f,'defaultTextFontName','Times New Roman');
                set(f,'defaultTextFontSize',8);
                
                set(f,'defaultLegendFontName','Times New Roman');
                set(f,'defaultLegendFontSize',8);


                
                % Define subplot positions manually as fractions of the figure height, shifted down by sgtitleSpace
                sgtitleSpace = 0.08;
                subplot1_position = [0.15, 0.84 - sgtitleSpace, 0.8, 0.15];
                subplot2_position = [0.15, 0.65 - sgtitleSpace, 0.8, 0.15];
                subplot3_position = [0.15, 0.45 - sgtitleSpace, 0.8, 0.15];
                subplot4_position = [0.15, 0.3 - sgtitleSpace, 0.8, 0.1];
                subplot5_position = [0.15, 0.0, 0.8, 0.25];

                
                ax1 = subplot('Position', subplot1_position);
                grid on 
                hold on
                box on
                plot(timeSimulink,gfmPQSimulink)
                plot(timeDMSIM,gfmPQDMSIM,'--')
                ylabel('$p,\,q$ (p.u.)',Interpreter='latex')
                xlim([timeDMSIM(1) timeDMSIM(end)])
                %ylim([-2 2])
                 xticklabels(ax1, {})
                sgtitle('Grid-forming Converter','FontSize',10,'FontWeight','bold')
                ax2 = subplot('Position', subplot2_position);
                grid on 
                hold on
                box on
                plot(timeSimulink,voltagePCCAmplitudeSimulink)
                plot(timeDMSIM,voltageAmplitudeDMSIM,'--')
                ylabel('$||\mathbf{v}_{dq}||$ (p.u.)',Interpreter='latex')
                xlim([timeDMSIM(1) timeDMSIM(end)])
               % ylim([-0.2 1.3])
                xticklabels(ax2, {})

                ax3 = subplot('Position', subplot3_position);
                grid on 
                hold on
                box on
                plot(timeSimulink,gfmCurrentAmplitudeSimulink)
                plot(timeDMSIM,gfmCurrentAmplitudeDMSIM,'--')
                ylabel('$||\mathbf{i}_{dq}||$ (p.u.)',Interpreter='latex')
                xlim([timeDMSIM(1) timeDMSIM(end)])
                %ylim([-0.2 2])
                xticklabels(ax3, {})


                ax4 = subplot('Position', subplot4_position);
                grid on 
                hold on
                box on
                plot(timeSimulink,frequencySimulink)
                plot(timeDMSIM,gfmFrequencyDMSIM,'--')
                ylabel('$f$ (Hz)',Interpreter='latex')
                %ylim([49.5 50.5])
                xlim([timeDMSIM(1) timeDMSIM(end)])
                xticklabels(ax4, {})

                ax5 = subplot('Position', subplot5_position);
                grid on 
                hold on
                box on
                plot(timeSimulink,deltaAngleSimulink)
                plot(timeDMSIM,deltaAngleDMSIM,'--')
                ylabel('$\delta (\circ)$',Interpreter='latex')
                xlim([timeDMSIM(1) timeDMSIM(end)])
                %ylim([1 2])

               % Set xlabel for the last plot
                xlabel(ax5, '$t$ (s)','FontSize',8,Interpreter='latex')
                legend('Nonlinear','Multilinear','Location','southoutside','Orientation','horizontal','NumColumns',2)
                
    %% grid-following
                figure%(2); clf;
                f = gcf;
                % Reserve space for sgtitle and increase figure height
                sgtitleSpace = 0.08;
                heightIEEEptWithTitle = heightIEEEpt / (1 - sgtitleSpace);
                set(f,'Units', 'points', 'Position', [100, 100, widthIEEEpt, heightIEEEptWithTitle],'Name','GFL Simulink MTI');

                set(f,'defaultLineLineWidth',1.2);
                
                set(f,'defaultAxesFontName','Times New Roman');
                set(f,'defaultAxesFontSize',8);
                
                set(f,'defaultTextFontName','Times New Roman');
                set(f,'defaultTextFontSize',8);
                
                set(f,'defaultLegendFontName','Times New Roman');
                set(f,'defaultLegendFontSize',8);


                
                % Define subplot positions manually as fractions of the figure height, shifted down by sgtitleSpace
                sgtitleSpace = 0.08;
                subplot1_position = [0.15, 0.84 - sgtitleSpace, 0.8, 0.15];
                subplot2_position = [0.15, 0.65 - sgtitleSpace, 0.8, 0.15];
                subplot3_position = [0.15, 0.45 - sgtitleSpace, 0.8, 0.15];
                subplot4_position = [0.15, 0.3 - sgtitleSpace, 0.8, 0.1];
                subplot5_position = [0.15, 0.0, 0.8, 0.25];

                sgtitle('Grid-following Converter','FontSize',10,'FontWeight','bold')
                ax1 = subplot('Position', subplot1_position);
                grid on 
                hold on
                box on
                plot(timeSimulink,gflPQSimulink)
                plot(timeDMSIM,gflPQDMSIM,'--')
                ylabel('$p,\,q$ (p.u.)',Interpreter='latex')
                xlim([timeDMSIM(1) timeDMSIM(end)])
                %ylim([-2 2])
                 xticklabels(ax1, {})

                ax2 = subplot('Position', subplot2_position);
                grid on 
                hold on
                box on
                plot(timeSimulink,voltagePCCAmplitudeSimulink)
                plot(timeDMSIM,voltageAmplitudeDMSIM,'--')
                ylabel('$||\mathbf{v}_{dq}||$ (p.u.)',Interpreter='latex')
                xlim([timeDMSIM(1) timeDMSIM(end)])
                %ylim([-0.2 1.3])
                xticklabels(ax2, {})

                ax3 = subplot('Position', subplot3_position);
                grid on 
                hold on
                box on
                plot(timeSimulink,gflCurrentAmplitudeSimulink)
                plot(timeDMSIM,gflCurrentAmplitudeDMSIM,'--')
                ylabel('$||\mathbf{i}_{dq}||$ (p.u.)',Interpreter='latex')
                xlim([timeDMSIM(1) timeDMSIM(end)])
                %ylim([-0.2 2])
                xticklabels(ax3, {})


                ax4 = subplot('Position', subplot4_position);
                grid on 
                hold on
                box on
                plot(timeSimulink,gflFrequencySimulink)
                plot(timeDMSIM,gflFrequencyDMSIM,'--')
                ylabel('$f$ (Hz)',Interpreter='latex')
                %ylim([49.5 50.5])
                xlim([timeDMSIM(1) timeDMSIM(end)])
                xticklabels(ax4, {})

                ax5 = subplot('Position', subplot5_position);
                grid on 
                hold on
                box on
                plot(timeSimulink,gfldeltaAngleSimulink)
                plot(timeDMSIM,gfldeltaAngleDMSIM,'--')
                ylabel('$\delta (\circ)$',Interpreter='latex')
                xlim([timeDMSIM(1) timeDMSIM(end)])
                %ylim([1 2])

               % Set xlabel for the last plot
                xlabel(ax5, '$t$ (s)','FontSize',8,Interpreter='latex')
                legend('Nonlinear','Multilinear','Location','southoutside','Orientation','horizontal','NumColumns',2)
                
               


%% 
               %  figure%(3); clf;
               %  f = gcf;
               %  set(f,'Units', 'points', 'Position', [100, 100,widthIEEEpt ,  0.8*heightIEEEpt],'Name','GFM GFL MTI');
               % 
               %  set(f,'defaultLineLineWidth',1.2);
               % 
               %  set(f,'defaultAxesFontName','Times New Roman');
               %  set(f,'defaultAxesFontSize',8);
               % 
               %  set(f,'defaultTextFontName','Times New Roman');
               %  set(f,'defaultTextFontSize',8);
               % 
               %  set(f,'defaultLegendFontName','Times New Roman');
               %  set(f,'defaultLegendFontSize',8);
               % 
               % 
               % 
               %  % Define subplot positions manually as fractions of the figure height
               %  subplot1_position = [0.15, 0.8, 0.8, 0.15];
               %  subplot2_position = [0.15, 0.7, 0.8, 0.15]; % Reduced height for the second plot
               %  subplot3_position = [0.15, 0.6, 0.8, 0.15];  
               %  subplot4_position =[0.15, 0.4, 0.8, 0.15] ; 
               %  subplot5_position = [0.15, 0, 0.8, 0.325];%[0.15, 0.0, 0.8, 0.25];  
               % 
               % 
               %  ax1 = subplot('Position', subplot1_position);
               %  grid on 
               %  hold on
               %  box on
               %  plot(timeDMSIM,gflPQDMSIM,'-')
               %  plot(timeDMSIM,gfmPQDMSIM,'--')
               %  ylabel('$P,\,Q$ (p.u.)',Interpreter='latex')
               %  xlim([timeDMSIM(1) timeDMSIM(end)])
               %  ylim([-2 2])
               %   xticklabels(ax1, {})
               %  ax3 = subplot('Position', subplot3_position);
               %  grid on 
               %  hold on
               %  box on
               % 
               %  plot(timeDMSIM,gflCurrentAmplitudeDMSIM,'-')
               %  plot(timeDMSIM,gfmCurrentAmplitudeDMSIM,'--')
               %  ylabel('$||\mathbf{I}_{dq}||$ (p.u.)',Interpreter='latex')
               %  xlim([timeDMSIM(1) timeDMSIM(end)])
               %  ylim([-0.2 2])
               %  xticklabels(ax3, {})
               % 
               % 
               %  ax4 = subplot('Position', subplot4_position);
               %  grid on 
               %  hold on
               %  box on
               %  plot(timeDMSIM,gflFrequencyDMSIM,'-')
               %  plot(timeDMSIM,gfmFrequencyDMSIM,'--')
               % 
               %  ylabel('$f$ (Hz)',Interpreter='latex')
               %  ylim([49.5 50.5])
               %  xlim([timeDMSIM(1) timeDMSIM(end)])
               %  xticklabels(ax4, {})
               % 
               %  ax5 = subplot('Position', subplot5_position);
               %  grid on 
               %  hold on
               %  box on
               %  plot(timeDMSIM,gfldeltaAngleDMSIM,'-')
               %  plot(timeDMSIM,deltaAngleDMSIM,'--')
               % 
               %  ylabel('$\delta (\circ)$',Interpreter='latex')
               %  xlim([timeDMSIM(1) timeDMSIM(end)])
               % 
               % % Set xlabel for the last plot
               %  xlabel(ax5, '$t$ (s)','FontSize',8,Interpreter='latex')
               %  legend('Grid-following converter','Grid-forming converter','Location','southoutside','Orientation','horizontal','NumColumns',2)
                
                

        
       