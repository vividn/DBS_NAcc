clear all; close all; hold all;
soma = importdata('soma.dat','\t',2); soma = soma.data;
lastnode = importdata('axon.dat','\t',2); lastnode = lastnode.data;

[~,somaPeakI] = findpeaks(soma(:,2),'MINPEAKHEIGHT',-20);
[~,nodePeakI] = findpeaks(lastnode(:,2),'MINPEAKHEIGHT',-20);

somaPeaks = soma(somaPeakI,1);
nodePeaks = lastnode(nodePeakI,1);


gold = [200 165 40] ./ 255.0;
maroon = [122 0 25] ./ 255.0;

plot(soma(:,1),soma(:,2),'LineWidth',2,'Color',gold);
legend('soma')
ylim([-90 10])
xlim([0 2000])
hline(-55,'--')
hline(-64,'--')

figure;
plot(lastnode(:,1),lastnode(:,2),'LineWidth',2,'Color',maroon);
ylim([-90 10])
xlim([0 2000])

ylabel('Membrane Potential (mV)');
xlabel('time (ms)');

legend('axon tip');

