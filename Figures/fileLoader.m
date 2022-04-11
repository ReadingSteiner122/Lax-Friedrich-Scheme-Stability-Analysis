
files  = dir('C:\Users\sanath\Documents\Projects\NPDE\Assignment-2\matlabplots\*.dat');

for i=1:length(files)
    load(files(i).name, '-ascii');
end

for i=1:length(files)
    figure(i);
    plot(files(i).name(:,1), files(i).name(:, 2));
    hold on;
    plot(Exacta2var501(:, 1), Exacta2var501(:, 2));
    hold off;
end
    
    