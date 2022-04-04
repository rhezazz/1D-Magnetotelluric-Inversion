%Program pemodelan inversi kurva sounding MT 1-D dengan
%menggunakan algoritma FA (Firefly Algorithm)
%Mohammad Rheza Zamani
tic
clear all;
clc;
%Data sintetik
R = [100 10 1000];
thk = [500 1000];
freq = logspace(-3,3,50);
T = 1./freq;
[app_sin, phase_sin] = modelMT(R, thk ,T);

%Definisi ruang model 
npop = 100; %Jumlah dari model  
nlayer = 3; %Jumlah lapisan 
nitr = 500; %Jumlah iterasi 
%Ruang pencarian diatur sebesar 5 kali dari model data  sintetik
%Batas bawah pencarian nilai resistivitas
LBR = [1 1 1];
%Batas atas pencarian nilai resistivitas
UBR = [500 50 5000];
%Batas bawah pencarian nilai ketebalan
LBT = [1 1];
%Batas atas pencarian nilai resistivitas
UBT = [2500 5000];
alpha = 0.2;
betha0 = 1;
gamma = 0.8;
damp = 0.99;
%Membuat model awal
%Membuat model awal acak
for ipop = 1 : npop
    rho(ipop , :) = LBR + rand*(UBR - LBR);
    thick(ipop, :) = LBT + rand*(UBT - LBT);
end

for ipop=1:npop
    [apparentResistivity, phase_baru]=modelMT(rho(ipop,:),thick(ipop,:),T);
     app_mod(ipop,:)=apparentResistivity;
     phase_mod(ipop,:)=phase_baru;
     
    [misfit]=misfitMT(app_sin,phase_sin,app_mod(ipop,:),phase_mod(ipop,:));
    E(ipop)=misfit;
end
%Inversi
for itr = 1 : nitr
    for i =  1 : npop
        j = randi(npop,1);
        while i == j
            j = randi(npop,1);
        end
        if E(i)<E(j)
            %Calculated Distance
            dr = norm((rho(i,:)-rho(j,:)));
            dt = norm((thick(i,:)-thick(j,:)));
            %Calculated new model with determine a new position
            %Random vector position for resistivity model
            for n1 = 1 : nlayer
                rho_baru(1,n1) = rho(i,n1) +betha0.*exp(-gamma*(dr)^2).*(rho(j,n1)-rho(i,n1))+ (alpha*(rand-0.5)*abs((UBR(n1)-LBR(n1))));
                if rho_baru(1,n1) < LBR(n1);
                     rho_baru(1,n1) = LBR(n1);
                end
                if rho_baru(1,n1) > UBR(n1);
                     rho_baru(1,n1) = UBR(n1);
                end
            end
            %Random vector position for thick model
            for n2 = 1 : (nlayer-1)
                thk_baru(1,n2) = thick(i,n2) +betha0.*exp(-gamma*(dt)^2).*(thick(j,n2)-thick(i,n2))+ (alpha*(rand-0.5)*abs((UBT(n2)-LBT(n2))));
                if thk_baru(1,n2) < LBT(n2);
                     thk_baru(1,n2) = LBT(n2);
                end
                if thk_baru(1,n2) > UBT(n2);
                    thk_baru(1,n2) = UBT(n2);
                end
            end 
        else
            rho_baru(1,:) = rho(i,:);
            thk_baru(1,:) = thick(i,:);
        end
        [apparentResistivity_baru, phase_baru]=modelMT(rho_baru,thk_baru,T);
        [err] = misfitMT(app_sin,phase_sin,apparentResistivity_baru, phase_baru);
        %Update model dan error jika lebih baik
        if err<E(i)
            rho(i,:) = rho_baru(1,:);
            thick(i,:) = thk_baru(1,:);
            app_mod(i,:) = apparentResistivity_baru(1,:);
            phase_mod(i,:) = phase_baru(1,:);
            E(i) = err;
        end
      
    end
     Emin = 100;
        for ipop = 1 : npop
        if E(ipop)< Emin
            Emin = E(ipop);
            rho_model = rho(ipop,:);
            thk_model = thick(ipop,:);
            app_model = app_mod(ipop,:);
            phase_model = phase_mod(ipop,:);
        end
    end
    Egen(itr)=Emin;
    alpha = alpha*damp;
end
time = toc
    %Persiapan Ploting
    rho_plot = [0 R];
    thk_plot = [0 cumsum(thk) max(thk)*10000];
    rhomod_plot = [0 rho_model];
    thkmod_plot = [0 cumsum(thk_model) max(thk_model)*10000];
    %Plotting
    figure(1)
    subplot(2, 2, 1)
    loglog(T,app_sin,'.b',T,app_model,'r','MarkerSize',12,'LineWidth',1.5);
    axis([10^-3 10^3 1 10^3]);
    legend({'Synthetic Data','Calculated Data'},'EdgeColor','none','Color','none','FontWeight','Bold','location','southeast');
    xlabel('Periods (s)','FontSize',12,'FontWeight','Bold');
    ylabel('App. Resistivity (Ohm.m)','FontSize',12,'FontWeight','Bold');
    title(['\bf \fontsize{10}\fontname{Times}Period (s) vs Apparent Resistivity (ohm.m)  || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)]);
    grid on
    
    subplot(2, 2, 3)
    loglog(T,phase_sin,'.b',T,phase_model,'r','MarkerSize',12,'LineWidth',1.5);
    axis([10^-3 10^3 0 90]);
    set(gca, 'YScale', 'linear');
    legend({'Synthetic Data','Calculated Data'},'EdgeColor','none','Color','none','FontWeight','Bold');
    xlabel('Periods (s)','FontSize',12,'FontWeight','Bold');
    ylabel('Phase (deg)','FontSize',12,'FontWeight','Bold');
    title(['\bf \fontsize{10}\fontname{Times}Period (s) vs Phase (deg)  || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)]);
    grid on
    
    subplot(2, 2, [2 4])
    stairs(rho_plot,thk_plot,'--r','Linewidth',1.5);
    hold on
    stairs(rhomod_plot ,thkmod_plot,'-b','Linewidth',2);
    hold off
     legend({'Synthetic Model','Calculated Model'},'EdgeColor','none','Color','none','FontWeight','Bold','Location','SouthWest');
    axis([1 10^4 0 5000]);
    xlabel('Resistivity (Ohm.m)','FontSize',12,'FontWeight','Bold');
    ylabel('Depth (m)','FontSize',12,'FontWeight','Bold');
    title(['\bf \fontsize{10}\fontname{Times}Model']);
    subtitle(['\rho_{1} = ',num2str(rho_model(1)),' || \rho_{2} = ',num2str(rho_model(2)),' || \rho_{3} = ',num2str(rho_model(3)),' || thick_{1} = ',num2str(thk_model(1)),' || thick_{2} = ',num2str(thk_model(2))],'FontWeight','bold')
    set(gca,'YDir','Reverse');
    set(gca, 'XScale', 'log');
    set(gcf, 'Position', get(0, 'Screensize'));
    grid on

%plot misfit
figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RSME','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
grid on