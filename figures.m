clc
clear
close all
load("data\results_article.mat") %5GNR 1024
clear label labe2

itnr="{\it N_r}";
itnp="{\itN_p}";
itdp="{\it d_p}";
itn="{\it N}";
itdr="{\it d_r}";

lineStyle=["-", "--",":"];
marker=["o","^"];
color=["#0072BD", "#D95319", "#7E2F8E"];
geoName=["c-URA","d-URA","d-ELSA"];

%% FIGURE 3 of the article
% figure 3(a)
figure
a=1; k=1; i=1;
plot((nR(k)^2)*(nP(k,:).^2),squeeze(r_mat(a,k,i,:)),'-*',Color=color(a))
hold on
label(i,k,a)=geoName(a);
for a=2:nG
    for k=1:length(nR)
        for i=1:length(dP)
            plot((nR(k)^2)*(nP(k,:).^2),squeeze(r_mat(a,k,i,:)),lineStyle(k)+marker(i),Color=color(a))
            label(i,k,a)=geoName(a)+","+itnr+"="+string(nR(k)^2)+","+itdp+"="+sprintf("%0.1f",dP(i)*lambda)+" m";
        end
    end
end
hold off
ylabel('Coverage radius (km)')

grid on
box on
ylim([0,25])
xlim([0,1200])

label=reshape(label,[1,nG*length(nR)*length(dP)]);
label=rmmissing(label);
xlabel('Total number of radiating elements ({\itN})');
legend(label,Location="eastoutside")
%
% figure 3(b)
figure
ylabel('Throughput (Mb/s)')
ylim([0,300])
yticks(0:75:300)
hold on
a=1; k=1; i=1;
plot((nR(k)^2)*(nP(k,:).^2),squeeze(mbps1_mat(a,k,i,:)),'-*',Color=color(a))
for a=2:nG
    for k=1:length(nR)
        for i=1:length(dP)
            plot((nR(k)^2)*(nP(k,:).^2),squeeze(mbps1_mat(a,k,i,:)),lineStyle(k)+marker(i),Color=color(a))
        end
    end
end
yyaxis right
a=1; k=1; i=1;
ylabel('SNR (dB)')
ylim([-10,30])
yticks(-10:10:30)
xlim([0,1200])
plot((nR(k)^2)*(nP(k,:).^2),squeeze(SNR(a,k,i,:)),'-x',Color='#196f3d')
hold off
% title('(b)','FontWeight','Normal')
grid on
box on

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = '#196f3d';
xlabel('Total number of radiating elements ({\itN})')
label(10)='SNR';
legend(label,Location="eastoutside")

% figure 3(c)
figure
a=1; k=1; i=1;
semilogy((nR(k)^2)*(nP(k,:).^2),squeeze(thpDensity1_mat(a,k,i,:)),'-*',Color=color(a))
hold on
for a=2:nG
    for k=1:length(nR)
        for i=1:length(dP)
            semilogy((nR(k)^2)*(nP(k,:).^2),squeeze(thpDensity1_mat(a,k,i,:)),lineStyle(k)+marker(i),Color=color(a))
        end
    end
end
hold off

ylabel('Throughput area density $\left(\frac{\mathrm{Mb/s}}{\mathrm{km^2}}\right)$','Interpreter','latex')
ylim([1e-2,1e3])
yticks([1e-2,1e-1,1e0,1e1,1e2,1e3])
xlim([0,1200])

grid on
ax = gca;
ytickformat('%.0f'); % Set format to decimals
ax.YAxis.Exponent = 0; % Remove exponent
box on

xlabel('Total number of radiating elements ({\itN})')
label(10)='SNR';
legend(label(1:9),Location="eastoutside")
%% FIGURE 5 of the article

% figure 5(b)
figure
geoName2=["c-URA (Noise only)","d-URA","d-ELSA"];
a=1; k=1; i=1;
semilogy((nR(k)^2)*(nP(k,:).^2),squeeze(thpDensity1_mat(a,k,i,:)),'-*',Color=color(a))
label2(i,k,a)=geoName2(a);
hold on
for a=3:3
    for k=1:length(nR)
        for i=1:length(dP)
            semilogy((nR(k)^2)*(nP(k,:).^2),squeeze(thpDensity2_mat(a,k,i,:)),lineStyle(k)+marker(i),Color=color(a))
            % label2(i,k,a)=sprintf("%s,%s=%d,%s=%0.1f m \nSINR= %1d dB",geoName2(a),itnr,nR(k)^2,itdp,dP(i)*lambda,SINR_lim);
            label2(i,k,a)=sprintf("%s,%s=%d,%s=%0.1f m",geoName2(a),itnr,nR(k)^2,itdp,dP(i)*lambda);
        end
    end
end
hold off
label2=reshape(label2,[1,nG*length(nR)*length(dP)]);
label2=rmmissing(label2);


ylabel('Throughput area density $\left(\frac{\mathrm{Mb/s}}{\mathrm{km^2}}\right)$','Interpreter','latex')
ylim([1e-2,1e1])
xlim([0,1200])
grid on
ax = gca;
% Set format to decimals
ytickformat('%.0f');
% Remove exponent
ax.YAxis.Exponent = 0;
box on
xlabel('Total number of radiating elements ({\itN})');
legend(label2,Location="eastoutside")

% figure 5(a)
figure
for k=1:length(nR)
    hold on
    for i=1:length(dP)
        plot((nR(k)^2)*(nP(k,:).^2),squeeze(iB_mat(k,i,:)),lineStyle(k)+marker(i),Color=color(a))
    end
end
hold off
grid on
ax = gca;
ytickformat('%.0f');
ax.YAxis.Exponent = 0;
box on

ylabel('Number of interfering beams')
xlim([0,1200])
ylim([0,700])
yticks(0:100:700)
xlabel('Total number of radiating elements ({\itN})')
legend(label2(2:end),Location="eastoutside")