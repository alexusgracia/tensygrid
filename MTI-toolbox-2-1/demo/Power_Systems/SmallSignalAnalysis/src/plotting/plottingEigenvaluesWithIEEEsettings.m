% 
% Plotting the figures of the paper according to IEEE standard
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
                figure(6)

                % Define the figure with specific width and height in points
                f=figure('Units', 'points', 'Position', [100, 100,widthIEEEpt ,  0.5*widthIEEEpt]);

                set(f,'defaultLineLineWidth',1.2);
                
                set(f,'defaultAxesFontName','Times New Roman');
                set(f,'defaultAxesFontSize',8);
                
                set(f,'defaultTextFontName','Times New Roman');
                set(f,'defaultTextFontSize',8);
                
                set(f,'defaultLegendFontName','Times New Roman');
                set(f,'defaultLegendFontSize',8);

                
                hold on
                plot(eigRef,'o')
                plot(eigenvaluesDSS,'x')
                grid on
                box on
                ylabel('Imaginary$(\lambda)$',Interpreter='latex')
                xlabel('Real$(\lambda)$','FontSize',8,Interpreter='latex')
                xlim([-40 5])
                
                legend({'Eigenvalues', 'Generalized Eigenvalues'}, 'Location', 'southoutside')


