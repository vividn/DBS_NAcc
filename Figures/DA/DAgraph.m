DA = importdata('DA.dat');
x = -90:-55;
plot(x,DA)
legend(num2str((1.0:0.1:1.4)'));
xlabel('membrane voltage (mV)')
ylabel('potassium and calcium ion conductivity (µS/cm²)')