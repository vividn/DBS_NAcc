clear all; close all; hold all;
soma = importdata('soma.dat','\t',2); soma = soma.data;
lastnode = importdata('lastnode.dat','\t',2); lastnode = lastnode.data;

gold = [255 204 51] ./ 255.0;
maroon = [122 0 25] ./ 255.0;

plot(soma(:,1),soma(:,2),'LineWidth',2,'Color',gold);
hold on;
plot(lastnode(:,1),lastnode(:,2),'LineWidth',2,'Color',maroon);

ylabel('Membrane Potential (V)');
xlabel('time (ms)');

legend('soma','axon tip');