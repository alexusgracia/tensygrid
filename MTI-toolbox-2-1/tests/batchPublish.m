files = dir(fullfile(fileparts(which("mssDoc.m")), '*.m'));
for i = 1:length(files)
    clearvars -EXCEPT files i;
    publish(files(i).name,'outputDir',strcat(files(i).folder,'\output'))
end