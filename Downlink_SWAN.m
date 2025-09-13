%% Downlink SWAN
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
Monte_Carlo = 1e4;
L = 1;
Side_Length = [10:10:200];
PASS1 = zeros(1,length(Side_Length));
PASS2 = zeros(1,length(Side_Length));
PASS3 = zeros(1,length(Side_Length));
PASS4 = zeros(1,length(Side_Length));
PASS5 = zeros(1,length(Side_Length));
PASS6 = zeros(1,length(Side_Length));
PASS7 = zeros(1,length(Side_Length));
PASS8 = zeros(1,length(Side_Length));
PASS9 = zeros(1,length(Side_Length));
PASS10 = zeros(1,length(Side_Length));
PASS11 = zeros(1,length(Side_Length));
PASS12 = zeros(1,length(Side_Length));
PASS13 = zeros(1,length(Side_Length));
PASS14 = zeros(1,length(Side_Length));
PASS15 = zeros(1,length(Side_Length));
parfor Monte = [1:1:Monte_Carlo]
    uy = rand(1)*20 - 10;
    ux0 = rand(1);
    cy = h^2 + uy^2;
    PASS1_tmp = zeros(1,length(Side_Length));
    PASS2_tmp = zeros(1,length(Side_Length));
    PASS3_tmp = zeros(1,length(Side_Length));
    PASS4_tmp = zeros(1,length(Side_Length));
    PASS5_tmp = zeros(1,length(Side_Length));
    PASS6_tmp = zeros(1,length(Side_Length));
    PASS7_tmp = zeros(1,length(Side_Length));
    PASS8_tmp = zeros(1,length(Side_Length));
    PASS9_tmp = zeros(1,length(Side_Length));
    PASS10_tmp = zeros(1,length(Side_Length));
    PASS11_tmp = zeros(1,length(Side_Length));
    PASS12_tmp = zeros(1,length(Side_Length));
    PASS13_tmp = zeros(1,length(Side_Length));
    PASS14_tmp = zeros(1,length(Side_Length));
    PASS15_tmp = zeros(1,length(Side_Length));
    for index = [1:1:length(Side_Length)]
        [index,Monte]
        Dx = Side_Length(index);
        ux = ux0 * Dx;
        %% SWAN - SS
        L = 1;
        M = Dx / L;
        if ux / L == M
            M0 = M;
        elseif mod(ux, L) == 0
            M0 = ux / L + 1;
        else
            M0 = ceil(ux / L);
        end
        d_previous = ux;
        distance_right = [];        
        while (d_previous + lambda / 2) <= (M0 * L)
            psi0 = (M0 - 1) * L;
            psi1 = d_previous + lambda / 2;
            d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
            d = d0 + mod(sqrt(cy) + neff * (ux - ((M0 - 1) * L)) - d0, lambda);
            Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
            v = ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1) - psi1;
            d_previous = psi1 + v;
            distance_right = [distance_right, d_previous];
        end
        d_previous = ux;
        distance_left = [];
        while (d_previous - lambda / 2) >= ((M0 - 1) * L)
            psi0 = (M0 - 1) * L;
            psi1 = d_previous - lambda / 2;
            d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
            d = d0 - mod(d0 - (sqrt(cy) + neff * (ux - ((M0 - 1) * L))), lambda);
            Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
            v = psi1 - ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1);
            d_previous = psi1 - v;
            distance_left = [distance_left, d_previous];
        end
        distance = [flip(distance_left), ux, distance_right];
        N = length(distance);
        PASS4_tmp(index) = log2(1 + power * eta / N / noise * abs(sum(1./sqrt(cy + (distance - ux).^2)...
            .*exp(-1*sqrt(-1)*2*pi/lambda*(sqrt((distance - ux).^2 + cy) + neff * (distance - (M0 - 1) * L)))))^2);
        PASS5_tmp(index) = log2(1 + power * eta / N / noise * abs(sum(1./sqrt(cy + (distance - ux).^2)))^2);
        PASS6_tmp(index) = log2(1 + power * eta / N / noise * abs(sum(1./sqrt(cy + (distance - ux).^2).*exp(-1*alpha*abs(distance - (M0 - 1) * L))))^2);
        %% Conventional PASS
        M0 = 1;
        L = Dx / M0;
        d_previous = ux;
        distance_right = ones(1, N - 1);        
        for n_index = [1:1:(N - 1)]
            psi0 = 0;
            psi1 = d_previous + lambda / 2;
            d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
            d = d0 + mod(sqrt(cy) + neff * (ux - ((M0 - 1) * L)) - d0, lambda);
            Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
            v = ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1) - psi1;
            d_previous = psi1 + v;
            distance_right(n_index) = d_previous;
        end
        d_previous = ux;
        distance_left = ones(1, N - 1);        
        for n_index = [1:1:(N - 1)]
            psi0 = 0;
            psi1 = d_previous - lambda / 2;
            d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
            d = d0 - mod(d0 - (sqrt(cy) + neff * (ux - ((M0 - 1) * L))), lambda);
            Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
            v = psi1 - ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1);
            d_previous = psi1 - v;
            distance_left(n_index) = d_previous;
        end
        N_right = N - 1 - length(find(distance_right > M0 * L));
        N_left = N - 1 - length(find(distance_left < ((M0 - 1) * L)));
        if N_right < (N - 1)/2
            N_left = N - 1 - N_right;
        elseif N_left < (N - 1)/2
            N_right = N - 1 - N_left;
        elseif (N_right >= (N - 1)/2)&&(N_left > (N - 1)/2)
            N_right = floor((N - 1)/2);
            N_left = N - 1 - N_right;
        elseif (N_right > (N - 1)/2)&&(N_left >= (N - 1)/2)
            N_left = floor((N - 1)/2);
            N_right = N - 1 - N_left;
        end
        distance = [flip(distance_left([1:1:N_left])), ux, distance_right([1:1:N_right])];
        PASS1_tmp(index) = log2(1 + power * eta / N / noise * abs(sum(1./sqrt(cy + (distance - ux).^2)...
            .*exp(-1*sqrt(-1)*2*pi/lambda*(sqrt((distance - ux).^2 + cy) + neff * (distance - (M0 - 1) * L)))))^2);
        PASS2_tmp(index) = log2(1 + power * eta / N / noise * abs(sum(1./sqrt(cy + (distance - ux).^2)))^2);
        PASS3_tmp(index) = log2(1 + power * eta / N / noise * abs(sum(1./sqrt(cy + (distance - ux).^2).*exp(-1*alpha*abs(distance - (M0 - 1) * L))))^2);        
        %% SWAN - SA
        L = 1;
        M = Dx / L;
        if ux / L == M
            M0 = M;
        elseif mod(ux, L) == 0
            M0 = ux / L + 1;
        else
            M0 = ceil(ux / L);
        end
        location = cell(M,1);
        d_previous = ux;
        distance_right = [];        
        while (d_previous + lambda / 2) <= (M0 * L)
            psi0 = (M0 - 1) * L;
            psi1 = d_previous + lambda / 2;
            d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
            d = d0 + mod(sqrt(cy) + neff * (ux - ((M0 - 1) * L)) - d0, lambda);
            Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
            v = ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1) - psi1;
            d_previous = psi1 + v;
            distance_right = [distance_right, d_previous];
        end
        d_previous = ux;
        distance_left = [];
        while (d_previous - lambda / 2) >= ((M0 - 1) * L)
            psi0 = (M0 - 1) * L;
            psi1 = d_previous - lambda / 2;
            d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
            d = d0 - mod(d0 - (sqrt(cy) + neff * (ux - ((M0 - 1) * L))), lambda);
            Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
            v = psi1 - ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1);
            d_previous = psi1 - v;
            distance_left = [distance_left, d_previous];
        end
        distance = [flip(distance_left), ux, distance_right];
        location{M0,1} = distance;
        for m_index = [(M0 + 1):1:M]
            d_previous = location{m_index - 1,1};
            d_previous = d_previous(end);
            distance_right = [];        
            while (d_previous + lambda / 2) <= (m_index * L)
                psi0 = (m_index - 1) * L;
                psi1 = d_previous + lambda / 2;
                d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
                d = d0 + mod(sqrt(cy) + neff * (ux - ((M0 - 1) * L)) - d0, lambda);
                Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
                v = ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1) - psi1;
                d_previous = psi1 + v;
                distance_right = [distance_right, d_previous];
            end    
            location{m_index,1} = distance_right;
        end
        for m_index = [(M0 - 1):(-1):1]
            d_previous = location{m_index + 1,1};
            d_previous = d_previous(1);
            distance_left = [];
            while (d_previous - lambda / 2) >= ((m_index - 1) * L)
                psi0 = (m_index - 1) * L;
                psi1 = d_previous - lambda / 2;
                d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
                d = d0 - mod(d0 - (sqrt(cy) + neff * (ux - ((M0 - 1) * L))), lambda);
                Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
                v = psi1 - ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1);
                d_previous = psi1 - v;
                distance_left = [distance_left, d_previous];
            end    
            location{m_index,1} = flip(distance_left);
        end
        PASS10_tmp0 = 0;
        PASS11_tmp0 = 0;
        PASS12_tmp0 = 0;
        N_total = 0;
        for m_index = [1:1:M]
            distance = location{m_index,1};
            N_total = N_total + length(distance);
            PASS10_tmp0 = PASS10_tmp0 + sum(1/sqrt(length(distance))./sqrt(cy + (distance - ux).^2)...
            .*exp(-1*sqrt(-1)*2*pi/lambda*(sqrt((distance - ux).^2 + cy) + neff * (distance - (m_index - 1) * L))));
            PASS11_tmp0 = PASS11_tmp0 + sum(1/sqrt(length(distance))./sqrt(cy + (distance - ux).^2));
            PASS12_tmp0 = PASS12_tmp0 + sum(1/sqrt(length(distance))./sqrt(cy + (distance - ux).^2).*exp(-1*alpha*abs(distance - (m_index - 1) * L))); 
        end
        PASS10_tmp(index) = log2(1 + power * eta / M / noise * abs(PASS10_tmp0)^2);
        PASS11_tmp(index) = log2(1 + power * eta / M / noise * abs(PASS11_tmp0)^2);
        PASS12_tmp(index) = log2(1 + power * eta / M / noise * abs(PASS12_tmp0)^2);
        %% Conventional PASS
        M0 = 1;
        L = Dx / M0;
%         N = N_total;
%         d_previous = ux;
%         distance_right = ones(1, N - 1);        
%         for n_index = [1:1:(N - 1)]
%             psi0 = 0;
%             psi1 = d_previous + lambda / 2;
%             d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
%             d = d0 + mod(sqrt(cy) + neff * (ux - ((M0 - 1) * L)) - d0, lambda);
%             Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
%             v = ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1) - psi1;
%             d_previous = psi1 + v;
%             distance_right(n_index) = d_previous;
%         end
%         d_previous = ux;
%         distance_left = ones(1, N - 1);        
%         for n_index = [1:1:(N - 1)]
%             psi0 = 0;
%             psi1 = d_previous - lambda / 2;
%             d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
%             d = d0 - mod(d0 - (sqrt(cy) + neff * (ux - ((M0 - 1) * L))), lambda);
%             Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
%             v = psi1 - ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1);
%             d_previous = psi1 - v;
%             distance_left(n_index) = d_previous;
%         end
%         N_right = N - 1 - length(find(distance_right > M0 * L));
%         N_left = N - 1 - length(find(distance_left < ((M
%         .0 - 1) * L)));
%         if N_right < (N - 1)/2
%             N_left = N - 1 - N_right;
%         elseif N_left < (N - 1)/2
%             N_right = N - 1 - N_left;
%         end
%         distance = [flip(distance_left([1:1:N_left])), ux, distance_right([1:1:N_right])];
        
        d_previous = ux;
        distance_right = [];        
        while (d_previous + lambda / 2) <= (M0 * L)
            psi0 = (M0 - 1) * L;
            psi1 = d_previous + lambda / 2;
            d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
            d = d0 + mod(sqrt(cy) + neff * (ux - ((M0 - 1) * L)) - d0, lambda);
            Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
            v = ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1) - psi1;
            d_previous = psi1 + v;
            distance_right = [distance_right, d_previous];
        end
        d_previous = ux;
        distance_left = [];
        while (d_previous - lambda / 2) >= ((M0 - 1) * L)
            psi0 = (M0 - 1) * L;
            psi1 = d_previous - lambda / 2;
            d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
            d = d0 - mod(d0 - (sqrt(cy) + neff * (ux - ((M0 - 1) * L))), lambda);
            Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
            v = psi1 - ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1);
            d_previous = psi1 - v;
            distance_left = [distance_left, d_previous];
        end
        distance = [flip(distance_left), ux, distance_right];
        N = length(distance);    
        PASS7_tmp(index) = log2(1 + power * eta / N / noise * abs(sum(1./sqrt(cy + (distance - ux).^2)...
            .*exp(-1*sqrt(-1)*2*pi/lambda*(sqrt((distance - ux).^2 + cy) + neff * (distance - (M0 - 1) * L)))))^2);
        PASS8_tmp(index) = log2(1 + power * eta / N / noise * abs(sum(1./sqrt(cy + (distance - ux).^2)))^2);
        PASS9_tmp(index) = log2(1 + power * eta / N / noise * abs(sum(1./sqrt(cy + (distance - ux).^2).*exp(-1*alpha*abs(distance - (M0 - 1) * L))))^2);        
        %% SWAN - SM
        L = 1;
        M = Dx / L;
        if ux / L == M
            M0 = M;
        elseif mod(ux, L) == 0
            M0 = ux / L + 1;
        else
            M0 = ceil(ux / L);
        end
        location = cell(M,1);
        d_previous = ux;
        distance_right = [];        
        while (d_previous + lambda / 2) <= (M0 * L)
            psi0 = (M0 - 1) * L;
            psi1 = d_previous + lambda / 2;
            d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
            d = d0 + mod(sqrt(cy) + neff * (ux - ((M0 - 1) * L)) - d0, lambda);
            Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
            v = ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1) - psi1;
            d_previous = psi1 + v;
            distance_right = [distance_right, d_previous];
        end
        d_previous = ux;
        distance_left = [];
        while (d_previous - lambda / 2) >= ((M0 - 1) * L)
            psi0 = (M0 - 1) * L;
            psi1 = d_previous - lambda / 2;
            d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
            d = d0 - mod(d0 - (sqrt(cy) + neff * (ux - ((M0 - 1) * L))), lambda);
            Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
            v = psi1 - ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1);
            d_previous = psi1 - v;
            distance_left = [distance_left, d_previous];
        end
        distance = [flip(distance_left), ux, distance_right];
        location{M0,1} = distance;
        for m_index = [(M0 + 1):1:M]
            d_previous0 = location{m_index - 1,1};
            d_previous = d_previous0(end);
            distance_right = [];        
            while (d_previous + lambda / 2) <= (m_index * L)
                psi0 = (m_index - 1) * L;
                psi1 = d_previous + lambda / 2;
                d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
                d = d0 + mod(sqrt(cy + (d_previous0(end) + lambda / 2 - ux)^2) + neff * ((d_previous0(end) + lambda / 2) - ((m_index - 1) * L)) - d0, lambda);
                Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
                v = ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1) - psi1;
                d_previous = psi1 + v;
                distance_right = [distance_right, d_previous];
            end    
            location{m_index,1} = distance_right;
        end
        for m_index = [(M0 - 1):(-1):1]
            d_previous0 = location{m_index + 1,1};
            d_previous = d_previous0(1);
            distance_left = [];
            while (d_previous - lambda / 2) >= ((m_index - 1) * L)
                psi0 = (m_index - 1) * L;
                psi1 = d_previous - lambda / 2;
                d0 = sqrt((psi1 - ux)^2 + cy) + neff * (psi1 - psi0);
                d = d0 - mod(d0 - (sqrt(cy + (d_previous0(end) - lambda / 2 - ux)^2) + neff * ((d_previous0(end) - lambda / 2) - ((m_index - 1) * L))), lambda);
                Delta = (ux - psi0)^2 * (neff^2) - 2 * d * neff * (ux - psi0) + cy * (neff^2 - 1) + d^2;
                v = psi1 - ((neff^2) * psi0 + d * neff - sqrt(Delta) - ux)/(neff^2 - 1);
                d_previous = psi1 - v;
                distance_left = [distance_left, d_previous];
            end    
            location{m_index,1} = flip(distance_left);
        end
        PASS13_tmp0 = 0;
        PASS14_tmp0 = 0;
        PASS15_tmp0 = 0;
        for m_index = [1:1:M]
            distance = location{m_index,1};
            PASS13_tmp0 = PASS13_tmp0 + abs(sum(1/sqrt(length(distance))./sqrt(cy + (distance - ux).^2)...
            .*exp(-1*sqrt(-1)*2*pi/lambda*(sqrt((distance - ux).^2 + cy) + neff * (distance - (m_index - 1) * L)))))^2;
            PASS14_tmp0 = PASS14_tmp0 + abs(sum(1/sqrt(length(distance))./sqrt(cy + (distance - ux).^2)))^2;
            PASS15_tmp0 = PASS15_tmp0 + abs(sum(1/sqrt(length(distance))./sqrt(cy + (distance - ux).^2).*exp(-1*alpha*abs(distance - (m_index - 1) * L))))^2; 
        end
        PASS13_tmp(index) = log2(1 + power * eta / noise * abs(PASS13_tmp0));
        PASS14_tmp(index) = log2(1 + power * eta / noise * abs(PASS14_tmp0));
        PASS15_tmp(index) = log2(1 + power * eta / noise * abs(PASS15_tmp0));
    end
    PASS1 = PASS1 + PASS1_tmp;
    PASS2 = PASS2 + PASS2_tmp;
    PASS3 = PASS3 + PASS3_tmp;
    PASS4 = PASS4 + PASS4_tmp;
    PASS5 = PASS5 + PASS5_tmp;
    PASS6 = PASS6 + PASS6_tmp;
    PASS7 = PASS7 + PASS7_tmp;
    PASS8 = PASS8 + PASS8_tmp;
    PASS9 = PASS9 + PASS9_tmp;
    PASS10 = PASS10 + PASS10_tmp;
    PASS11 = PASS11 + PASS11_tmp;
    PASS12 = PASS12 + PASS12_tmp;
    PASS13 = PASS13 + PASS13_tmp;
    PASS14 = PASS14 + PASS14_tmp;
    PASS15 = PASS15 + PASS15_tmp;
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
plot(Side_Length,PASS8/Monte_Carlo,'-v','MarkerSize',MarkerSize_n,'Color',[53,144,58]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS9/Monte_Carlo,'--v','MarkerSize',MarkerSize_n,'Color',[53,144,58]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS10/Monte_Carlo,'--d','MarkerSize',MarkerSize_n,'Color',[53,144,58]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS11/Monte_Carlo,'-d','MarkerSize',MarkerSize_n,'Color',[53,144,58]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS12/Monte_Carlo,'--o','MarkerSize',MarkerSize_n,'Color',[53,144,58]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS13/Monte_Carlo,'--s','MarkerSize',MarkerSize_n,'Color',[53,144,58]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS14/Monte_Carlo,'--x','MarkerSize',MarkerSize_n,'Color',[53,144,58]/255,'LineWidth',LineWidth_n);
hold on;
plot(Side_Length,PASS15/Monte_Carlo,'--p','MarkerSize',MarkerSize_n,'Color',[53,144,58]/255,'LineWidth',LineWidth_n);
hold on;
xlim([1,201]);
Data = [PASS1;PASS2;PASS3;PASS4;PASS5;PASS6;PASS7;PASS8;PASS9;PASS10;PASS11;PASS12;PASS13;PASS14;PASS15];
save Downlink_SWAN_Average_Rate_Side_Length Data
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