%% Uplink SWAN
h = 3; % height
f = 28*1e9; % carrier frequency
c = 299792458; % speed of light
lambda = c/f; % wavelength
k0 = 2*pi/lambda; % wavenumber
neff = 1.4;
eta = (lambda/(4*pi))^2;
power = 10;
noise = 10^(-90/10); % noise power
kappa = 0.08;
alpha = log(10) * kappa / 20;
Monte_Carlo = 1000000;
L = 1;
Side_Length = [1:5:101];
PASS1 = zeros(1,length(Side_Length));
PASS2 = zeros(1,length(Side_Length));
PASS3 = zeros(1,length(Side_Length));
PASS4 = zeros(1,length(Side_Length));
PASS5 = zeros(1,length(Side_Length));
PASS6 = zeros(1,length(Side_Length));
PASS7 = zeros(1,length(Side_Length));
PASS8 = zeros(1,length(Side_Length));
PASS9 = zeros(1,length(Side_Length));
for index = [1:1:length(Side_Length)]
    Dx = Side_Length(index)*L;
    M = Dx / L;
    for Monte = [1:1:Monte_Carlo]
        uy = rand(1)*20 - 10;
        ux = Side_Length(index)*rand(1);
        cy = h^2 + uy^2;
        %% Conventional PASS
        PASS1(index) = PASS1(index) + log2(1 + (sqrt(eta) / sqrt(cy))^2 * power / noise);
        PASS2(index) = PASS2(index) + log2(1 + (sqrt(eta) / sqrt(cy) * exp(-1*alpha*abs(ux - 0))).^2 * power / noise);
        %% SWAN - SS
        PASS3(index) = PASS3(index) + log2(1 + (sqrt(eta) / sqrt(cy))^2 * power / noise);
        ux1 = mod(ux,L);
        PASS4(index) = PASS4(index) + log2(1 + (sqrt(eta) / sqrt(cy) * exp(-1*alpha*abs(ux1 - 0))).^2 * power / noise);
        %% SWAN - SA
        if ux / L == M
            M0 = M;
        elseif mod(ux, L) == 0
            M0 = ux / L + 1;
        else
            M0 = ceil(ux / L);
        end
        d_previous = ux;
        M_right = [(M0 + 1):1:M];
        distance_right = ones(1, length(M_right));
        for m_index = [1:1:length(M_right)]
            m = M_right(m_index);
            psi0 = (m - 1) * L;
            psi1 = max([d_previous + lambda / 2, psi0]);
            d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
            d = d0 + mod(sqrt(cy) + neff * (ux - ((M0 - 1) * L)) - d0, lambda);
            Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
            v = ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1) - psi1;
            d_previous = psi1 + v;
            distance_right(m_index) = d_previous;
        end
        d_previous = ux;
        M_left = [(M0 - 1):(-1):1];
        distance_left = ones(1, length(M_left));
        for m_index = [1:1:length(M_left)]
            m = M_left(m_index);
            psi0 = (m - 1) * L;
            psi1 = min([d_previous - lambda / 2, psi0 + L]);
            d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
            d = d0 - mod(d0 - (sqrt(cy) + neff * (ux - ((M0 - 1) * L))), lambda);
            Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
            v = psi1 - ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1);
            d_previous = psi1 - v;
            distance_left(m_index) = d_previous;
        end
        distance = [flip(distance_left), ux, distance_right];
        PASS5(index) = PASS5(index) + log2(1 + power * eta / M / noise * abs(sum(1./sqrt(cy + (distance - ux).^2)...
            .*exp(-1*sqrt(-1)*2*pi/lambda*(sqrt((distance - ux).^2 + cy) + neff * (distance - ([1:1:M] - 1)*L)))))^2);
        PASS6(index) = PASS6(index) + log2(1 + power * eta / M / noise * abs(sum(1./sqrt(cy + (distance - ux).^2)))^2);
        PASS7(index) = PASS7(index) + log2(1 + power * eta / M / noise * abs(sum(1./sqrt(cy + (distance - ux).^2).*exp(-1*alpha*abs(distance - ([1:1:M] - 1)*L))))^2);
        %% SWAN - SM
        if ux / L == M
            M0 = M;
        elseif mod(ux, L) == 0
            M0 = ux / L + 1;
        else
            M0 = ceil(ux / L);
        end
        d_previous = ux;
        M_right = [(M0 + 1):1:M];
        distance_right = ones(1, length(M_right));
        for m_index = [1:1:length(M_right)]
            m = M_right(m_index);
            psi0 = (m - 1) * L;
            psi1 = max([d_previous + lambda / 2, psi0]);
            v = 0;
            d_previous = psi1 + v;
            distance_right(m_index) = d_previous;
        end
        d_previous = ux;
        M_left = [(M0 - 1):(-1):1];
        distance_left = ones(1, length(M_left));
        for m_index = [1:1:length(M_left)]
            m = M_left(m_index);
            psi0 = (m - 1) * L;
            psi1 = min([d_previous - lambda / 2, psi0 + L]);
            v = 0;
            d_previous = psi1 - v;
            distance_left(m_index) = d_previous;
        end
        distance = [flip(distance_left), ux, distance_right];
        PASS8(index) = PASS8(index) + log2(1 + power * eta / noise * abs(sum(1./(cy + (distance - ux).^2))));
        PASS9(index) = PASS9(index) + log2(1 + power * eta / noise * abs(sum((1./sqrt(cy + (distance - ux).^2).*exp(-1*alpha*abs(distance - ([1:1:M] - 1)*L))).^2)));
    end
end
MarkerSize_n = 14;
LineWidth_n = 3;
plot(Side_Length,PASS1/Monte_Carlo,'-o','MarkerSize',MarkerSize_n,'Color',[53,144,58]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS2/Monte_Carlo,'-s','MarkerSize',MarkerSize_n,'Color',[53,144,58]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS3/Monte_Carlo,'-','MarkerSize',MarkerSize_n,'Color',[190,42,44]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS4/Monte_Carlo,'--','MarkerSize',MarkerSize_n,'Color',[32,110,158]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS5/Monte_Carlo,'-x','MarkerSize',MarkerSize_n,'Color',[32,110,158]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS6/Monte_Carlo,'-+','MarkerSize',MarkerSize_n,'Color',[53,144,58]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS7/Monte_Carlo,'--+','MarkerSize',MarkerSize_n,'Color',[53,144,58]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS8/Monte_Carlo,'-+','MarkerSize',MarkerSize_n,'Color',[53,144,58]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS9/Monte_Carlo,'--+','MarkerSize',MarkerSize_n,'Color',[53,144,58]/255,'LineWidth',LineWidth_n);
hold on;
xlim([1,101]);
Data = [PASS1;PASS2;PASS3;PASS4;PASS5;PASS6;PASS7;PASS8;PASS9];
save Uplink_SWAN_Average_Rate_Side_Length Data
% xlim([1,401]);
% Data = [Data;PASS1;PASS2;PASS3;PASS4;PASS5];
% ylabel('Uplink Received SNR','Fontsize',25, 'Fontname', 'Times New Roman');
% xlabel('Number of Segments','Fontsize',25,'Fontname', 'Times New Roman');
% h = legend('Conventional PASS','SWAN');
% set(h,'location','northeast','Fontsize',22,'Fontname', 'Times New Roman','ItemTokenSize',[81,54],'box','off');
% set(gca,'FontSize',25, 'Fontname', 'Times New Roman');
% set(gca,'XMinorTick','on');
% set(gca,'YMinorTick','on');
% set(gca,'LineWidth',1.5);
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 12 6];
% print('Uplink_Received_SNR_SA_Number_Segment.eps','-depsc','-r600');