% 
% Plotting the figures of the paper according to IEEE standard
% KauChr
% 01.05.2025
%

clf
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


%%                
               
                set(0,'defaultLineLineWidth',1.2);

%%
                figure(8)%;clf;
                f=gcf;
                
                % Define the figure with specific width and height in points
                set(f,'Units', 'points', 'Position', [100, 100,widthIEEEpt ,  0.5*widthIEEEpt]);

                set(f,'defaultLineLineWidth',1.2);
                
                set(f,'defaultAxesFontName','Times New Roman');
                set(f,'defaultAxesFontSize',8);
                
                set(f,'defaultTextFontName','Times New Roman');
                set(f,'defaultTextFontSize',8);
                
                set(f,'defaultLegendFontName','Times New Roman');
                set(f,'defaultLegendFontSize',8);

                
                hold on
              

                % Define a summer colormap
                colorMap=summer(length(trajectoryGeneralizedEigenvalues(1,:)));
                
                % Preallocate a cell array for legend entries
                legendEntries = cell(length(trajectoryGeneralizedEigenvalues(1,:)), 1);
                
             
                hold on
                for k=1:length(trajectoryGeneralizedEigenvalues(1,:))
                    
                   plot(real(trajectoryGeneralizedEigenvalues(:,k)),imag(trajectoryGeneralizedEigenvalues(:,k)), 'o', ...
                        'MarkerFaceColor', colorMap(k, :), ...
                        'MarkerEdgeColor', 'k', 'MarkerSize', 8)
                             
                    plot(real(trajectoryReferenceEigenvalues(:,k)),imag(trajectoryReferenceEigenvalues(:,k)), 'x', 'MarkerSize', 8)
                
                end
                

                grid on
             
               % Add legend for marker meanings
               legend({'Generalized Eigenvalues', 'Eigenvalues'}, 'Location', 'southoutside')
             
                box on
                ylabel('Imaginary$(\lambda)$',Interpreter='latex')
                xlabel('Real$(\lambda)$','FontSize',8,Interpreter='latex')
                xlim([-40 5])
                 
            

%%

