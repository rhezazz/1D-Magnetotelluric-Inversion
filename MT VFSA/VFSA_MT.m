%Program pemodelan inversi kurva sounding MT 1-D dengan
%menggunakan algoritma Very Fast Simulated Annealing
%Mohammad Rheza Zamani
clear all;
clc;
%Data sintetik
R = [100 10 1000];
thk = [500 1500];
freq = logspace(-3,3,50);
T = 1./freq;
[app_sin, phase_sin] = modelMT(R, thk ,T);

%Definisi ruang model 
nlayer = 3; %Jumlah lapisan 
nitr = 200; %Jumlah iterasi 
%Ruang pencarian diatur sebesar 5 kali dari model data  sintetik
%Batas bawah pencarian nilai resistivitas
LBR = [1 1 1];
%Batas atas pencarian nilai resistivitas
UBR = [200 20 2000];
%Batas bawah pencarian nilai ketebalan
LBT = [1 1];
%Batas atas pencarian nilai resistivitas
UBT = [1000 3000];
Temp = 1;
dec = 1;
%Membuat model awal acak
rho1(1 , :) = LBR + rand*(UBR - LBR);
thick1(1, :) = LBT + rand*(UBT - LBT);
%Menghitung misfit, apparent resistivity, dan phase model awal
[apparentResistivity1, phase1]=modelMT(rho1(1,:),thick1(1,:),T);
app_mod1(1,:)=apparentResistivity1;
phase_mod1(1,:)=phase1;
     
[misfit1]=misfitMT(app_sin,phase_sin,app_mod1(1,:),phase_mod1(1,:));
E1=misfit1;
for itr = 1 : nitr
    rho_int(1 , :) = LBR + rand*(UBR - LBR);
    thick_int(1, :) = LBT + rand*(UBT - LBT);
    ui = rand;
    yi = sign(ui-0.5)*Temp*((((1 + (1/Temp)))^abs(2*ui-1))-1);
    rho2(1 , :) = rho_int + yi*(UBR - LBR);
    thick2(1, :) = thick_int + yi*(UBT - LBT);
    [apparentResistivity2, phase2]=modelMT(rho2(1,:),thick2(1,:),T);
    app_mod2(1,:)=apparentResistivity2;
    phase_mod2(1,:)=phase2;
    [misfit2]=misfitMT(app_sin,phase_sin,app_mod2(1,:),phase_mod2(1,:));
    E2=misfit2;
    delta_E = E2 -E1;
    
    if delta_E < 0
        rho1 = rho2;
        thick1 = thick2;
         E1 = E2;
    else
        P = exp((-delta_E)/Temp);
        if  P >= rand
           rho1 = rho2;
           thick1 = thick2;
           E1 = E2;
        end
    end
    [apparentResistivity_new, phase_new]=modelMT(rho1(1,:),thick1(1,:),T);
    Egen(itr)=E1;
    Temp = Temp*exp(-dec*(itr)^(1/(2*nlayer)-1));
    Temperature(itr) = Temp;

    %Persiapan Ploting
    rho_plot = [0 R];
    thk_plot = [0 cumsum(thk) max(thk)*10000];
    rhomod_plot = [0 rho1];
    thkmod_plot = [0 cumsum(thick1) max(thick1)*10000];
    %Plotting
    figure(1)
    subplot(2, 2, 1)
    loglog(T,app_sin,'.b',T,apparentResistivity_new,'r','MarkerSize',12,'LineWidth',1.5);
    axis([10^-3 10^3 1 10^3]);
    legend({'Synthetic Data','Calculated Data'},'EdgeColor','none','Color','none','FontWeight','Bold');
    xlabel('Periods (s)','FontSize',12,'FontWeight','Bold');
    ylabel('App. Resistivity (Ohm.m)','FontSize',12,'FontWeight','Bold');
    title(['\bf \fontsize{10}\fontname{Times}Period (s) vs Apparent Resistivity (ohm.m)  || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)]);
    grid on
    
    subplot(2, 2, 3)
    loglog(T,phase_sin,'.b',T,phase_new,'r','MarkerSize',12,'LineWidth',1.5);
    axis([10^-3 10^3 0 90]);
    set(gca, 'YScale', 'linear');
    legend({'Synthetic Data','Calculated Data'},'EdgeColor','none','Color','none','FontWeight','Bold');
    xlabel('Periods (s)','FontSize',12,'FontWeight','Bold');
    ylabel('Phase (deg)','FontSize',12,'FontWeight','Bold');
    title(['\bf \fontsize{10}\fontname{Times}Period (s) vs Phase (deg)  || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)]);
    grid on
    
    subplot(2, 2, [2 4])
    stairs(rho_plot,thk_plot,'--b','Linewidth',1.5);
    hold on
    stairs(rhomod_plot ,thkmod_plot,'-r','Linewidth',2);
    hold off
    legend({'Synthetic Model','Calculated Model'},'EdgeColor','none','Color','none','FontWeight','Bold','Location','SouthEast');
    axis([1 10^4 0 5000]);
    xlabel('Resistivity (Ohm.m)','FontSize',12,'FontWeight','Bold');
    ylabel('Depth (m)','FontSize',12,'FontWeight','Bold');
    title(['\bf \fontsize{10}\fontname{Times}Model']);
    subtitle(['\rho_{1} = ',num2str(rho1(1)),' || \rho_{2} = ',num2str(rho1(2)),' || \rho_{3} = ',num2str(rho1(3)),' || thick_{1} = ',num2str(thick1(1)),' || thick_{2} = ',num2str(thick1(2))],'FontWeight','bold')
    set(gca,'YDir','Reverse');
    set(gca, 'XScale', 'log');
    set(gcf, 'Position', get(0, 'Screensize'));
    grid on
end

%plot misfit
figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RSME','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
set(gcf, 'Position', get(0, 'Screensize'));
grid on

%Plot Temperature
figure(3)
plot(1:nitr,Temperature,'b','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('Temperature','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Penurunan Temperature ');
set(gcf, 'Position', get(0, 'Screensize'));
grid on