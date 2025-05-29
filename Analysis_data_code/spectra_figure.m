%% load data 
clc
close all
clear all
load spectra_data_trimmed

% for wave case W1, data include relative velocity (slip velocity) of 
% each particle when z_p'/h <-0.3. 
% for still water, particle velocity for z_p'/h<-0.3 are given


%% compute spectra using FFT 
clc
close all

fs = 20;

% N1
for kk=1:10
    clear xx
    xx = rltv_nylon_332_x(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_x_N1(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

for kk=1:10
    clear xx
    xx = rltv_nylon_332_z(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_z_N1(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end
freq_N1 = freq;

% N2 
for kk=1:10
    clear xx
    xx = rltv_nylon_18_x(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_x_N2(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

for kk=1:10
    clear xx
    xx = rltv_nylon_18_z(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_z_N2(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end
freq_N2 = freq;

% N3
for kk=1:10
    clear xx
    xx = rltv_nylon_14_x(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_x_N3(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

for kk=1:10
    clear xx
    xx = rltv_nylon_14_z(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_z_N3(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end
freq_N3 = freq;
% N4
for kk=1:10
    clear xx
    xx = rltv_nylon12_6_x(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_x_N4(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

for kk=1:10
    clear xx
    xx = rltv_nylon12_6_z(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_z_N4(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end
freq_N4 = freq;

% T1

for kk=1:10
    clear xx
    xx = rltv_torlon_18_x(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_x_T1(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

for kk=1:10
    clear xx
    xx = rltv_torlon_18_z(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_z_T1(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

freq_T1 = freq;

% still
%N1

for kk=1:10
    clear xx
    xx = still_nylon_332_x(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_x_N1_still(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

for kk=1:10
    clear xx
    xx = still_nylon_332_z(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_z_N1_still(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

% N2 
for kk=1:10
    clear xx
    xx = still_nylon_18_x(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_x_N2_still(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

for kk=1:10
    clear xx
    xx = still_nylon_18_z(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_z_N2_still(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

% N3
for kk=1:10
    clear xx
    xx = still_nylon_14_x(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_x_N3_still(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

for kk=1:10
    clear xx
    xx = still_nylon_14_z(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_z_N3_still(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

% N4
for kk=1:10
    clear xx
    xx = still_nylon12_6_x(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_x_N4_still(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

for kk=1:10
    clear xx
    xx = still_nylon12_6_z(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_z_N4_still(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

% T1

for kk=1:10
    clear xx
    xx = still_torlon_18_x(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_x_T1_still(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

for kk=1:10
    clear xx
    xx = still_torlon_18_z(:,kk);
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_fin_z_T1_still(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

%% bootstrapping
clc
close all
nboot =2000;

for kk=1:size(psd_fin_x_N1,1)
    m = bootstrp(nboot,@mean,psd_fin_x_N1(kk,:));
    psd_bst_x_N1_W1(1,kk) = mean(m);
    ci_x_N1_W1(:,kk) = bootci(nboot,@mean,psd_fin_x_N1(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_z_N1(kk,:));
    psd_bst_z_N1_W1(1,kk) = mean(m);
    ci_z_N1_W1(:,kk) = bootci(nboot,@mean,psd_fin_z_N1(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_x_N1_still(kk,:));
    psd_bst_x_N1_still(1,kk) = mean(m);
    ci_x_N1_still(:,kk) = bootci(nboot,@mean,psd_fin_x_N1_still(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_z_N1_still(kk,:));
    psd_bst_z_N1_still(1,kk) = mean(m);
    ci_z_N1_still(:,kk) = bootci(nboot,@mean,psd_fin_z_N1_still(kk,:))-mean(m);
end



for kk=1:size(psd_fin_x_N2,1)
    m = bootstrp(nboot,@mean,psd_fin_x_N2(kk,:));
    psd_bst_x_N2_W1(1,kk) = mean(m);
    ci_x_N2_W1(:,kk) = bootci(nboot,@mean,psd_fin_x_N2(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_z_N2(kk,:));
    psd_bst_z_N2_W1(1,kk) = mean(m);
    ci_z_N2_W1(:,kk) = bootci(nboot,@mean,psd_fin_z_N2(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_x_N2_still(kk,:));
    psd_bst_x_N2_still(1,kk) = mean(m);
    ci_x_N2_still(:,kk) = bootci(nboot,@mean,psd_fin_x_N2_still(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_z_N2_still(kk,:));
    psd_bst_z_N2_still(1,kk) = mean(m);
    ci_z_N2_still(:,kk) = bootci(nboot,@mean,psd_fin_z_N2_still(kk,:))-mean(m);
end



for kk=1:size(psd_fin_x_N3,1)
    m = bootstrp(nboot,@mean,psd_fin_x_N3(kk,:));
    psd_bst_x_N3_W1(1,kk) = mean(m);
    ci_x_N3_W1(:,kk) = bootci(nboot,@mean,psd_fin_x_N3(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_z_N3(kk,:));
    psd_bst_z_N3_W1(1,kk) = mean(m);
    ci_z_N3_W1(:,kk) = bootci(nboot,@mean,psd_fin_z_N3(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_x_N3_still(kk,:));
    psd_bst_x_N3_still(1,kk) = mean(m);
    ci_x_N3_still(:,kk) = bootci(nboot,@mean,psd_fin_x_N3_still(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_z_N3_still(kk,:));
    psd_bst_z_N3_still(1,kk) = mean(m);
    ci_z_N3_still(:,kk) = bootci(nboot,@mean,psd_fin_z_N3_still(kk,:))-mean(m);
end



for kk=1:size(psd_fin_x_N4,1)
    m = bootstrp(nboot,@mean,psd_fin_x_N4(kk,:));
    psd_bst_x_N4_W1(1,kk) = mean(m);
    ci_x_N4_W1(:,kk) = bootci(nboot,@mean,psd_fin_x_N4(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_z_N4(kk,:));
    psd_bst_z_N4_W1(1,kk) = mean(m);
    ci_z_N4_W1(:,kk) = bootci(nboot,@mean,psd_fin_z_N4(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_x_N4_still(kk,:));
    psd_bst_x_N4_still(1,kk) = mean(m);
    ci_x_N4_still(:,kk) = bootci(nboot,@mean,psd_fin_x_N4_still(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_z_N4_still(kk,:));
    psd_bst_z_N4_still(1,kk) = mean(m);
    ci_z_N4_still(:,kk) = bootci(nboot,@mean,psd_fin_z_N4_still(kk,:))-mean(m);
end



for kk=1:size(psd_fin_x_T1,1)
    m = bootstrp(nboot,@mean,psd_fin_x_T1(kk,:));
    psd_bst_x_T1_W1(1,kk) = mean(m);
    ci_x_T1_W1(:,kk) = bootci(nboot,@mean,psd_fin_x_T1(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_z_T1(kk,:));
    psd_bst_z_T1_W1(1,kk) = mean(m);
    ci_z_T1_W1(:,kk) = bootci(nboot,@mean,psd_fin_z_T1(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_x_T1_still(kk,:));
    psd_bst_x_T1_still(1,kk) = mean(m);
    ci_x_T1_still(:,kk) = bootci(nboot,@mean,psd_fin_x_T1_still(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,psd_fin_z_T1_still(kk,:));
    psd_bst_z_T1_still(1,kk) = mean(m);
    ci_z_T1_still(:,kk) = bootci(nboot,@mean,psd_fin_z_T1_still(kk,:))-mean(m);
end

%% Figure 11: Normalized ensemble mean spectra
clc
close all
lin_wid =0.8;
dia = [0.0024;0.0032;0.0064;0.006;0.0032];
vg = [0.0529;0.0611;0.0865;0.042;0.1116];
y_min = 0;
y_max = 0.81;
x_min =0;
x_max = 0.2;
delta = 0;

f_wave = 9.56/(2*pi);


figure(1)
h(1)=subplot(231);

errorbar(freq_N1*dia(1)/vg(1)-delta,psd_bst_x_N1_W1/(dia(1)*vg(1))/10^6,ci_x_N1_W1(1,:)/(dia(1)*vg(1))/10^6,...
    ci_x_N1_W1(1,:)/(dia(1)*vg(1))/10^6,'-r','linewidth',lin_wid);
hold on
errorbar(freq_N1*dia(1)/vg(1)-delta,psd_bst_z_N1_W1/(dia(1)*vg(1))/10^6,ci_z_N1_W1(1,:)/(dia(1)*vg(1))/10^6,...
    ci_z_N1_W1(1,:)/(dia(1)*vg(1))/10^6,'-g','linewidth',lin_wid);
errorbar(freq_N1*dia(1)/vg(1)+delta,psd_bst_x_N1_still/(dia(1)*vg(1))/10^6,ci_x_N1_still(1,:)/(dia(1)*vg(1))/10^6,...
    ci_x_N1_still(1,:)/(dia(1)*vg(1))/10^6,'-b','linewidth',lin_wid);
errorbar(freq_N1*dia(1)/vg(1)+delta,psd_bst_z_N1_still/(dia(1)*vg(1))/10^6,ci_z_N1_still(1,:)/(dia(1)*vg(1))/10^6,...
    ci_z_N1_still(1,:)/(dia(1)*vg(1))/10^6,'-k','linewidth',lin_wid);

ylim([y_min y_max])
% xlabel('$fd_p/v_g$','interpreter','latex');
ylabel('$P(f) /v^{\prime}_g d_p$','interpreter','latex');
xlim([x_min x_max])

xticks(freq_N1(2:2:end)*dia(1)/vg(1))
xtickformat('%.3f')
plot(ones(2,1)*f_wave*dia(1)/vg(1),[0 1],'--k','linewidth',lin_wid)
grid on

h(2)=subplot(232);

errorbar(freq_N2*dia(2)/vg(2)-delta,psd_bst_x_N2_W1/(dia(2)*vg(2))/10^6,ci_x_N2_W1(1,:)/(dia(2)*vg(2))/10^6,...
    ci_x_N2_W1(1,:)/(dia(2)*vg(2))/10^6,'-r','linewidth',lin_wid);
hold on
errorbar(freq_N2*dia(2)/vg(2)-delta,psd_bst_z_N2_W1/(dia(2)*vg(2))/10^6,ci_z_N2_W1(1,:)/(dia(2)*vg(2))/10^6,...
    ci_z_N2_W1(1,:)/(dia(2)*vg(2))/10^6,'-g','linewidth',lin_wid);
errorbar(freq_N2*dia(2)/vg(2)+delta,psd_bst_x_N2_still/(dia(2)*vg(2))/10^6,ci_x_N2_still(1,:)/(dia(2)*vg(2))/10^6,...
    ci_x_N2_still(1,:)/(dia(2)*vg(2))/10^6,'-b','linewidth',lin_wid);
errorbar(freq_N2*dia(2)/vg(2)+delta,psd_bst_z_N2_still/(dia(2)*vg(2))/10^6,ci_z_N2_still(1,:)/(dia(2)*vg(2))/10^6,...
    ci_z_N2_still(1,:)/(dia(2)*vg(2))/10^6,'-k','linewidth',lin_wid);

ylim([y_min y_max])
% xlabel('$fd_p/v_g$','interpreter','latex');
ylabel('$P(f) /v^{\prime}_g d_p$','interpreter','latex');
xlim([x_min x_max])
xticks(freq_N2(2:end)*dia(2)/vg(2))
xtickformat('%.3f')
grid on
plot(ones(2,1)*f_wave*dia(2)/vg(2),[0 1],'--k','linewidth',lin_wid)


h(3)=subplot(233);

errorbar(freq_N3*dia(3)/vg(3)-delta,psd_bst_x_N3_W1/(dia(3)*vg(3))/10^6,ci_x_N3_W1(1,:)/(dia(3)*vg(3))/10^6,...
    ci_x_N3_W1(1,:)/(dia(3)*vg(3))/10^6,'-r','linewidth',lin_wid);
hold on
errorbar(freq_N3*dia(3)/vg(3)-delta,psd_bst_z_N3_W1/(dia(3)*vg(3))/10^6,ci_z_N3_W1(1,:)/(dia(3)*vg(3))/10^6,...
    ci_z_N3_W1(1,:)/(dia(3)*vg(3))/10^6,'-g','linewidth',lin_wid);
errorbar(freq_N3*dia(3)/vg(3)+delta,psd_bst_x_N3_still/(dia(3)*vg(3))/10^6,ci_x_N3_still(1,:)/(dia(3)*vg(3))/10^6,...
    ci_x_N3_still(1,:)/(dia(3)*vg(3))/10^6,'-b','linewidth',lin_wid);
errorbar(freq_N3*dia(3)/vg(3)+delta,psd_bst_z_N3_still/(dia(3)*vg(3))/10^6,ci_z_N3_still(1,:)/(dia(3)*vg(3))/10^6,...
    ci_z_N3_still(1,:)/(dia(3)*vg(3))/10^6,'-k','linewidth',lin_wid);

ylim([y_min y_max])
% xlabel('$fd_p/v_g$','interpreter','latex');
ylabel('$P(f) /v^{\prime}_g d_p$','interpreter','latex');
xlim([x_min x_max])
xlabel('$fd_p/v^{\prime}_g$','interpreter','latex');

xticks(freq_N3(2:end)*dia(3)/vg(3))
xtickformat('%.3f')
grid on

plot(ones(2,1)*f_wave*dia(3)/vg(3),[0 1],'--k','linewidth',lin_wid)

h(4)=subplot(234);

errorbar(freq_N4*dia(4)/vg(4)-delta,psd_bst_x_N4_W1/(dia(4)*vg(4))/10^6,ci_x_N4_W1(1,:)/(dia(4)*vg(4))/10^6,...
    ci_x_N4_W1(1,:)/(dia(4)*vg(4))/10^6,'-r','linewidth',lin_wid);
hold on
errorbar(freq_N4*dia(4)/vg(4)-delta,psd_bst_z_N4_W1/(dia(4)*vg(4))/10^6,ci_z_N4_W1(1,:)/(dia(4)*vg(4))/10^6,...
    ci_z_N4_W1(1,:)/(dia(4)*vg(4))/10^6,'-g','linewidth',lin_wid);
errorbar(freq_N4*dia(4)/vg(4)+delta,psd_bst_x_N4_still/(dia(4)*vg(4))/10^6,ci_x_N4_still(1,:)/(dia(4)*vg(4))/10^6,...
    ci_x_N4_still(1,:)/(dia(4)*vg(4))/10^6,'-b','linewidth',lin_wid);
errorbar(freq_N4*dia(4)/vg(4)+delta,psd_bst_z_N4_still/(dia(4)*vg(4))/10^6,ci_z_N4_still(1,:)/(dia(4)*vg(4))/10^6,...
    ci_z_N4_still(1,:)/(dia(4)*vg(4))/10^6,'-k','linewidth',lin_wid);

ylim([y_min y_max])
% xlabel('$fd_p/v_g$','interpreter','latex');
xlabel('$fd_p/v^{\prime}_g$','interpreter','latex');
ylabel('$P(f) /v^{\prime}_g d_p$','interpreter','latex');
xlim([x_min x_max])

xticks(freq_N4(2:end)*dia(4)/vg(4))
xtickformat('%.3f')

grid on
plot(ones(2,1)*f_wave*dia(4)/vg(4),[0 1],'--k','linewidth',lin_wid)

h(5)=subplot(235);

errorbar(freq_T1*dia(5)/vg(5)-delta,psd_bst_x_T1_W1/(dia(5)*vg(5))/10^6,ci_x_T1_W1(1,:)/(dia(5)*vg(5))/10^6,...
    ci_x_T1_W1(1,:)/(dia(5)*vg(5))/10^6,'-r','linewidth',lin_wid);
hold on
errorbar(freq_T1*dia(5)/vg(5)-delta,psd_bst_z_T1_W1/(dia(5)*vg(5))/10^6,ci_z_T1_W1(1,:)/(dia(5)*vg(5))/10^6,...
    ci_z_T1_W1(1,:)/(dia(5)*vg(5))/10^6,'-g','linewidth',lin_wid);
errorbar(freq_T1*dia(5)/vg(5)+delta,psd_bst_x_T1_still/(dia(5)*vg(5))/10^6,ci_x_T1_still(1,:)/(dia(5)*vg(5))/10^6,...
    ci_x_T1_still(1,:)/(dia(5)*vg(5))/10^6,'-b','linewidth',lin_wid);
errorbar(freq_T1*dia(5)/vg(5)+delta,psd_bst_z_T1_still/(dia(5)*vg(5))/10^6,ci_z_T1_still(1,:)/(dia(5)*vg(5))/10^6,...
    ci_z_T1_still(1,:)/(dia(5)*vg(5))/10^6,'-k','linewidth',lin_wid);

ylim([y_min y_max])
% xlabel('$fd_p/v_g$','interpreter','latex');
xlabel('$fd_p/v^{\prime}_g$','interpreter','latex');
ylabel('$P(f) /v^{\prime}_g d_p$','interpreter','latex');
xlim([x_min x_max])
xticks(freq_T1(2:end)*dia(5)/vg(5))
xtickformat('%.3f')
grid on
plot(ones(2,1)*f_wave*dia(5)/vg(5),[0 1],'--k','linewidth',lin_wid)

set(findall(gcf,'-property','FontSize'),'FontSize',12)
set(gcf,'position',[100 100 1300 500])

text(h(1),0.004,0.71,'(a)','interpreter','latex','FontSize',14)
text(h(2),0.004,0.71,'(b)','interpreter','latex','FontSize',14)
text(h(3),0.004,0.71,'(c)','interpreter','latex','FontSize',14)
text(h(4),0.004,0.71,'(d)','interpreter','latex','FontSize',14)
text(h(5),0.004,0.71,'(e)','interpreter','latex','FontSize',14)
set(h(1),'position',[0.1 0.57 0.22 0.35])
set(h(2),'position',[0.4 0.57 0.22 0.35])
set(h(3),'position',[0.7 0.57 0.22 0.35])
set(h(4),'position',[0.1 0.13 0.22 0.35])
set(h(5),'position',[0.4 0.13 0.22 0.35])

legend(h(1),{'$(v{^\prime}-u{^\prime})_x^*$ (W1)','$(v{^\prime}-u{^\prime})_z^*$ (W1)',...
    '${v{^\prime}_x}^*$ (Q)','${v{^\prime}_z}^*$ (Q)'},...
    'interpreter','latex','position',[0.7 0.22 0.2 0.14],'Fontsize',14,'numcolumns',1)
% Optional save image
% saveas(figure(1),'spectra_wave_still.jpg')
% saveas(figure(1),'spectra_wave_still.eps','epsc')
% saveas(figure(1),'spectra_wave_still.fig')

%%
clc
close all
clear all

load spectra_N2_raw_data

% for all wave and quiescent condition relative velocity of N2 paricle
% below z'_p/h <-0.3 are provided

%% separate data into low-frequency oblique (T) and chaotic regime (2T)
clc
close all

still_2T = [1 3 5 9 11 12 15 16 19 20 21 22 23 24 25];
still_T =[2 4 6 7 8 10 13 14 17 18];

W1_2T = [1 3 4 5 6 7 10 14 17 19 20 22 23 24 25];
W1_T =[2 8 9 11 12 13 15 18 21];

W3_2T =[2 3 5 10 12 13 15 16 18 19 20 21 22 24];
W3_T = [1 4 6 7 8 9 11 14 17 23 25];

W2_2T =[1 2 3 6 8 9 12 15 16 18 19 21 23 25];
W2_T = [4 5 7 10 11 13 14 17 20 22 24];
%%
clc
close all
fs = 20;
dt =1/fs;
auto_leng =32;
gamma = 1.12;
prtc_dia = 25.4/8/1000;
g = 9.81;
vg = sqrt((gamma-1)*g*prtc_dia)
for kk=1:size(still_nylon_18_x,2)
    auto_x_still(kk,:) = autocorr(still_nylon_18_x(:,kk),auto_leng);
    auto_z_still(kk,:) = autocorr(still_nylon_18_z(:,kk),auto_leng);
    auto_x_W1(kk,:) = autocorr(W1_nylon_18_x(:,kk),auto_leng);
    auto_z_W1(kk,:) = autocorr(W1_nylon_18_z(:,kk),auto_leng);
    auto_x_W2(kk,:) = autocorr(W2_nylon_18_x(:,kk),auto_leng);
    auto_z_W2(kk,:) = autocorr(W2_nylon_18_z(:,kk),auto_leng);
    auto_x_W3(kk,:) = autocorr(W3_nylon_18_x(:,kk),auto_leng);
    auto_z_W3(kk,:) = autocorr(W3_nylon_18_z(:,kk),auto_leng);    
end

mean_auto_x_still = mean(auto_x_still,1);
mean_auto_z_still = mean(auto_z_still,1);
mean_auto_x_W1 = mean(auto_x_W1,1);
mean_auto_z_W1 = mean(auto_z_W1,1);
mean_auto_x_W2 = mean(auto_x_W2,1);
mean_auto_z_W2 = mean(auto_z_W2,1);
mean_auto_x_W3 = mean(auto_x_W3,1);
mean_auto_z_W3 = mean(auto_z_W3,1);

tt_auto = (0:auto_leng)/fs;

%% spectra
clc
close all

dia_prtc=0.0032; %(m)
rho_sg = 1.12;
vg=sqrt((rho_sg-1)*9.81*dia_prtc); %(m/s)
N2_vs = -94.6;
fs=20;
data_ini = 1;
data_end = 33;

% mode 1: mean subtrated, mode2: mean remain
mde = 1;

for kk=1:length(still_2T)
    clear xx
    test_kk = still_2T(kk);
    xx = still_nylon_18_x(data_ini:data_end,test_kk);
    if mde==1
        xx=xx-mean(xx); 
    else
        xx=xx;
    end
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    still_2T_psd_fin_x(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

fin_psd_x_still_2T = mean(still_2T_psd_fin_x,2);

for kk=1:length(still_T)
    clear xx
    test_kk = still_T(kk);
    xx = still_nylon_18_x(data_ini:data_end,test_kk);
    if mde==1
        xx=xx-mean(xx); 
    else
        xx=xx;
    end
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    still_T_psd_fin_x(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

fin_psd_x_still_T = mean(still_T_psd_fin_x,2);

for kk=1:length(still_2T)
    clear xx
    test_kk = still_2T(kk);
    xx = still_nylon_18_z(data_ini:data_end,test_kk);
    xx= xx-mean(xx);
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    still_2T_psd_fin_z(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end
%
fin_psd_z_still_2T = mean(still_2T_psd_fin_z,2);

for kk=1:length(still_T)
    clear xx
    test_kk = still_T(kk);
    xx = still_nylon_18_z(data_ini:data_end,test_kk);
    xx= xx-mean(xx);
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    still_T_psd_fin_z(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

fin_psd_z_still_T = mean(still_T_psd_fin_z,2);

%W1



for kk=1:length(W1_2T)
    clear xx
    test_kk = W1_2T(kk);
    xx = W1_nylon_18_x(data_ini:data_end,test_kk);
        if mde==1
        xx=xx-mean(xx); 
    else
        xx=xx;
    end
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    W1_2T_psd_fin_x(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end
%
fin_psd_x_W1_2T = mean(W1_2T_psd_fin_x,2);

for kk=1:length(W1_T)
    clear xx
    test_kk = W1_T(kk);
    xx = W1_nylon_18_x(data_ini:data_end,test_kk);
        if mde==1
        xx=xx-mean(xx); 
    else
        xx=xx;
    end 
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    W1_T_psd_fin_x(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

fin_psd_x_W1_T = mean(W1_T_psd_fin_x,2);

for kk=1:length(W1_2T)
    clear xx
    test_kk = W1_2T(kk);
    xx = W1_nylon_18_z(data_ini:data_end,test_kk);
    xx= xx-mean(xx);
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    W1_2T_psd_fin_z(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end
%
fin_psd_z_W1_2T = mean(W1_2T_psd_fin_z,2);

for kk=1:length(W1_T)
    clear xx
    test_kk = W1_T(kk);
    xx = W1_nylon_18_z(data_ini:data_end,test_kk);
    xx= xx-mean(xx);
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    W1_T_psd_fin_z(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

fin_psd_z_W1_T = mean(W1_T_psd_fin_z,2);

%W2



for kk=1:length(W2_2T)
    clear xx
    test_kk = W2_2T(kk);
    xx = W2_nylon_18_x(data_ini:data_end,test_kk);
    if mde==1
        xx=xx-mean(xx); 
    else
        xx=xx;
    end
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    W2_2T_psd_fin_x(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end
%
fin_psd_x_W2_2T = mean(W2_2T_psd_fin_x,2);

for kk=1:length(W2_T)
    clear xx
    test_kk = W2_T(kk);
    xx = W2_nylon_18_x(data_ini:data_end,test_kk);
    if mde==1
        xx=xx-mean(xx); 
    else
        xx=xx;
    end
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    W2_T_psd_fin_x(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

fin_psd_x_W2_T = mean(W2_T_psd_fin_x,2);

for kk=1:length(W2_2T)
    clear xx
    test_kk = W2_2T(kk);
    xx = W2_nylon_18_z(data_ini:data_end,test_kk);
    xx= xx-mean(xx);
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    W2_2T_psd_fin_z(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end
%
fin_psd_z_W2_2T = mean(W2_2T_psd_fin_z,2);

for kk=1:length(W2_T)
    clear xx
    test_kk = W2_T(kk);
    xx = W2_nylon_18_z(data_ini:data_end,test_kk);
    xx= xx-mean(xx);
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    W2_T_psd_fin_z(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

fin_psd_z_W2_T = mean(W2_T_psd_fin_z,2);

%W3


for kk=1:length(W3_2T)
    clear xx
    test_kk = W3_2T(kk);
    xx = W3_nylon_18_x(data_ini:data_end,test_kk);
    if mde==1
        xx=xx-mean(xx); 
    else
        xx=xx;
    end
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    W3_2T_psd_fin_x(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end
%
fin_psd_x_W3_2T = mean(W3_2T_psd_fin_x,2);

for kk=1:length(W3_T)
    clear xx
    test_kk = W3_T(kk);
    xx = W3_nylon_18_x(data_ini:data_end,test_kk);
    if mde==1
        xx=xx-mean(xx); 
    else
        xx=xx;
    end
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    W3_T_psd_fin_x(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

fin_psd_x_W3_T = mean(W3_T_psd_fin_x,2);

for kk=1:length(W3_2T)
    clear xx
    test_kk = W3_2T(kk);
    xx = W3_nylon_18_z(data_ini:data_end,test_kk);
    xx= xx-mean(xx);
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    W3_2T_psd_fin_z(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end
%
fin_psd_z_W3_2T = mean(W3_2T_psd_fin_z,2);

for kk=1:length(W3_T)
    clear xx
    test_kk = W3_T(kk);
    xx = W3_nylon_18_z(data_ini:data_end,test_kk);
    xx= xx-mean(xx);
    num_sam = length(xx);   
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    W3_T_psd_fin_z(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

fin_psd_z_W3_T = mean(W3_T_psd_fin_z,2);


%figure

clc
close all

figure(1)
plot(freq*dia_prtc/vg,fin_psd_x_W1_T/(dia_prtc*vg)/10^6,'-ro')
hold on

plot(freq*dia_prtc/vg,fin_psd_x_W2_T/(dia_prtc*vg)/10^6,'-gd')
plot(freq*dia_prtc/vg,fin_psd_x_W3_T/(dia_prtc*vg)/10^6,'-b>')
plot(freq*dia_prtc/vg,fin_psd_x_still_T/(dia_prtc*vg)/10^6,'-k<')

    
plot(freq*dia_prtc/vg,fin_psd_z_W1_T/(dia_prtc*vg)/10^6,'--ro')
hold on
plot(freq*dia_prtc/vg,fin_psd_z_W2_T/(dia_prtc*vg)/10^6,'--gd')
plot(freq*dia_prtc/vg,fin_psd_z_W3_T/(dia_prtc*vg)/10^6,'--b>')
plot(freq*dia_prtc/vg,fin_psd_z_still_T/(dia_prtc*vg)/10^6,'--k<')
xlim([0 0.2])
if mde==1
    ylim([0 0.4])
else
    ylim([0 1.2])
end
xlabel('$fd_p/v^{\prime}_g$','interpreter','latex');
ylabel('$P(f) /v^{\prime}_g d_p$','interpreter','latex');
legend('$(v^{\prime}-u^{\prime})^*_x$ (W1)','$(v^{\prime}-u^{\prime})^*_x$ (W2)','$(v^{\prime}-u^{\prime})^*_x$ (W3)',...
    '$v^{\prime*}_x$ (Q)',...
    '$(v^{\prime}-u^{\prime})^*_z$ (W1)','$(v^{\prime}-u^{\prime})^*_z$ (W2)','$(v^{\prime}-u^{\prime})^*_z$ (W3)','$v^{\prime*}_z$ (Q)',...
    'interpreter','latex','numcolumns',2)
set(findall(gcf,'-property','FontSize'),'FontSize',12)


figure(2)

plot(freq*dia_prtc/vg,fin_psd_x_W1_2T/(dia_prtc*vg)/10^6,'-ro')
hold on
plot(freq*dia_prtc/vg,fin_psd_x_W2_2T/(dia_prtc*vg)/10^6,'-gd')
plot(freq*dia_prtc/vg,fin_psd_x_W3_2T/(dia_prtc*vg)/10^6,'-b>')
plot(freq*dia_prtc/vg,fin_psd_x_still_2T/(dia_prtc*vg)/10^6,'-k<')
xlim([0 0.2])
    
plot(freq*dia_prtc/vg,fin_psd_z_W1_2T/(dia_prtc*vg)/10^6,'--ro')
hold on

plot(freq*dia_prtc/vg,fin_psd_z_W2_2T/(dia_prtc*vg)/10^6,'--gd')
plot(freq*dia_prtc/vg,fin_psd_z_W3_2T/(dia_prtc*vg)/10^6,'--b>')
plot(freq*dia_prtc/vg,fin_psd_z_still_2T/(dia_prtc*vg)/10^6,'--k<')
if mde==1
    ylim([0 0.4])
else
    ylim([0 1.2])
end
xlabel('$fd_p/v^{\prime}_g$','interpreter','latex');
ylabel('$P(f) /v^{\prime}_g d_p$','interpreter','latex');
legend('$(v^{\prime}-u^{\prime})^*_x$ (W1)','$(v^{\prime}-u^{\prime})^*_x$ (W2)','$(v^{\prime}-u^{\prime})^*_x$ (W3)',...
    '$v^{\prime*}_x$ (Q)',...
    '$(v^{\prime}-u^{\prime})^*_z$ (W1)','$(v^{\prime}-u^{\prime})^*_z$ (W2)','$(v^{\prime}-u^{\prime})^*_z$ (W3)','$v^{\prime*}_z$ (Q)',...
    'interpreter','latex','numcolumns',2)
set(findall(gcf,'-property','FontSize'),'FontSize',12)

%% bootstrapping 
clc
nboot = 2000;

%2T (chaotic)
for kk=1:length(freq)
    m = bootstrp(nboot,@mean,still_2T_psd_fin_x(kk,:));
    psd_bst_x_still_2T(1,kk) = mean(m);
    ci_x_still_2T(:,kk) = bootci(nboot,@mean,still_2T_psd_fin_x(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,still_2T_psd_fin_z(kk,:));
    psd_bst_z_still_2T(1,kk) = mean(m);
    ci_z_still_2T(:,kk) = bootci(nboot,@mean,still_2T_psd_fin_z(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,W1_2T_psd_fin_x(kk,:));
    psd_bst_x_W1_2T(1,kk) = mean(m);
    ci_x_W1_2T(:,kk) = bootci(nboot,@mean,W1_2T_psd_fin_x(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,W1_2T_psd_fin_z(kk,:));
    psd_bst_z_W1_2T(1,kk) = mean(m);
    ci_z_W1_2T(:,kk) = bootci(nboot,@mean,W1_2T_psd_fin_z(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,W2_2T_psd_fin_x(kk,:));
    psd_bst_x_W2_2T(1,kk) = mean(m);
    ci_x_W2_2T(:,kk) = bootci(nboot,@mean,W2_2T_psd_fin_x(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,W2_2T_psd_fin_z(kk,:));
    psd_bst_z_W2_2T(1,kk) = mean(m);
    ci_z_W2_2T(:,kk) = bootci(nboot,@mean,W2_2T_psd_fin_z(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,W3_2T_psd_fin_x(kk,:));
    psd_bst_x_W3_2T(1,kk) = mean(m);
    ci_x_W3_2T(:,kk) = bootci(nboot,@mean,W3_2T_psd_fin_x(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,W3_2T_psd_fin_z(kk,:));
    psd_bst_z_W3_2T(1,kk) = mean(m);
    ci_z_W3_2T(:,kk) = bootci(nboot,@mean,W3_2T_psd_fin_z(kk,:))-mean(m);
    
end

% T (oblique oscillating)

for kk=1:length(freq)
    m = bootstrp(nboot,@mean,still_T_psd_fin_x(kk,:));
    psd_bst_x_still_T(1,kk) = mean(m);
    ci_x_still_T(:,kk) = bootci(nboot,@mean,still_T_psd_fin_x(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,still_T_psd_fin_z(kk,:));
    psd_bst_z_still_T(1,kk) = mean(m);
    ci_z_still_T(:,kk) = bootci(nboot,@mean,still_T_psd_fin_z(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,W1_T_psd_fin_x(kk,:));
    psd_bst_x_W1_T(1,kk) = mean(m);
    ci_x_W1_T(:,kk) = bootci(nboot,@mean,W1_T_psd_fin_x(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,W1_T_psd_fin_z(kk,:));
    psd_bst_z_W1_T(1,kk) = mean(m);
    ci_z_W1_T(:,kk) = bootci(nboot,@mean,W1_T_psd_fin_z(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,W2_T_psd_fin_x(kk,:));
    psd_bst_x_W2_T(1,kk) = mean(m);
    ci_x_W2_T(:,kk) = bootci(nboot,@mean,W2_T_psd_fin_x(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,W2_T_psd_fin_z(kk,:));
    psd_bst_z_W2_T(1,kk) = mean(m);
    ci_z_W2_T(:,kk) = bootci(nboot,@mean,W2_T_psd_fin_z(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,W3_T_psd_fin_x(kk,:));
    psd_bst_x_W3_T(1,kk) = mean(m);
    ci_x_W3_T(:,kk) = bootci(nboot,@mean,W3_T_psd_fin_x(kk,:))-mean(m);
    
    m = bootstrp(nboot,@mean,W3_T_psd_fin_z(kk,:));
    psd_bst_z_W3_T(1,kk) = mean(m);
    ci_z_W3_T(:,kk) = bootci(nboot,@mean,W3_T_psd_fin_z(kk,:))-mean(m);
    
end
%% figure 14: spectra for Osillating oblique and chaotic regimes

clc
close all
lin_wid =0.5;
lin_wid2=1;
delta = 0;
dia_prtc = 0.0032;
vg = 0.0611;

h(1)=subplot(211);
errorbar(freq*dia_prtc/vg-1.5*delta,psd_bst_x_W1_T/(dia_prtc*vg)/10^6,ci_x_W1_T(1,:)/(dia_prtc*vg)/10^6,...
    ci_x_W1_T(2,:)/(dia_prtc*vg)/10^6,'-or','linewidth',lin_wid);
hold on
errorbar(freq*dia_prtc/vg-0.5*delta,psd_bst_x_W2_T/(dia_prtc*vg)/10^6,ci_x_W2_T(1,:)/(dia_prtc*vg)/10^6,...
    ci_x_W2_T(2,:)/(dia_prtc*vg)/10^6,'-dg','linewidth',lin_wid);
errorbar(freq*dia_prtc/vg+0.5*delta,psd_bst_x_W3_T/(dia_prtc*vg)/10^6,ci_x_W3_T(1,:)/(dia_prtc*vg)/10^6,...
    ci_x_W3_T(2,:)/(dia_prtc*vg)/10^6,'->b','linewidth',lin_wid);
errorbar(freq*dia_prtc/vg+1.5*delta,psd_bst_x_still_T/(dia_prtc*vg)/10^6,ci_x_still_T(1,:)/(dia_prtc*vg)/10^6,...
    ci_x_still_T(2,:)/(dia_prtc*vg)/10^6,'-<k','linewidth',lin_wid);


errorbar(freq*dia_prtc/vg-1.5*delta,psd_bst_z_W1_T/(dia_prtc*vg)/10^6,ci_z_W1_T(1,:)/(dia_prtc*vg)/10^6,...
    ci_z_W1_T(2,:)/(dia_prtc*vg)/10^6,':or','linewidth',lin_wid2);
errorbar(freq*dia_prtc/vg-0.5*delta,psd_bst_z_W2_T/(dia_prtc*vg)/10^6,ci_z_W2_T(1,:)/(dia_prtc*vg)/10^6,...
    ci_z_W2_T(2,:)/(dia_prtc*vg)/10^6,':dg','linewidth',lin_wid2);
errorbar(freq*dia_prtc/vg+0.5*delta,psd_bst_z_W3_T/(dia_prtc*vg)/10^6,ci_z_W3_T(1,:)/(dia_prtc*vg)/10^6,...
    ci_z_W3_T(2,:)/(dia_prtc*vg)/10^6,':>b','linewidth',lin_wid2);
errorbar(freq*dia_prtc/vg+1.5*delta,psd_bst_z_still_T/(dia_prtc*vg)/10^6,ci_z_still_T(1,:)/(dia_prtc*vg)/10^6,...
    ci_z_still_T(2,:)/(dia_prtc*vg)/10^6,':<k','linewidth',lin_wid2);

if mde==1
    ylim([0 .45])
    yticks(0:0.1:0.4);
else
    ylim([0 1.5])
end
xlabel('$fd_p/v^{\prime}_g$','interpreter','latex');
ylabel('$P(f) /v^{\prime}_g d_p$','interpreter','latex');
xlim([0 0.15])

f_wave = [9.56 5.49 7.56]/(2*pi);
plot(ones(2,1)*f_wave(1)*dia_prtc/vg,[0 0.5],'--r','linewidth',lin_wid)
plot(ones(2,1)*f_wave(2)*dia_prtc/vg,[0 0.5],'--g','linewidth',lin_wid)
plot(ones(2,1)*f_wave(3)*dia_prtc/vg,[0 0.5],'--b','linewidth',lin_wid)


xticks(freq*dia_prtc/vg)
grid on
xtickformat('%.3f')

h(2)=subplot(212);
errorbar(freq*dia_prtc/vg-1.5*delta,psd_bst_x_W1_2T/(dia_prtc*vg)/10^6,ci_x_W1_2T(1,:)/(dia_prtc*vg)/10^6,...
    ci_x_W1_2T(2,:)/(dia_prtc*vg)/10^6,'-or','linewidth',lin_wid);
hold on
errorbar(freq*dia_prtc/vg-1.5*delta,psd_bst_z_W1_2T/(dia_prtc*vg)/10^6,ci_z_W1_2T(1,:)/(dia_prtc*vg)/10^6,...
    ci_z_W1_2T(2,:)/(dia_prtc*vg)/10^6,':or','linewidth',lin_wid2);

errorbar(freq*dia_prtc/vg-0.5*delta,psd_bst_x_W2_2T/(dia_prtc*vg)/10^6,ci_x_W2_2T(1,:)/(dia_prtc*vg)/10^6,...
    ci_x_W2_2T(2,:)/(dia_prtc*vg)/10^6,'-dg','linewidth',lin_wid);
errorbar(freq*dia_prtc/vg-0.5*delta,psd_bst_z_W2_2T/(dia_prtc*vg)/10^6,ci_z_W2_2T(1,:)/(dia_prtc*vg)/10^6,...
    ci_z_W2_2T(2,:)/(dia_prtc*vg)/10^6,':dg','linewidth',lin_wid2);
errorbar(freq*dia_prtc/vg+0.5*delta,psd_bst_x_W3_2T/(dia_prtc*vg)/10^6,ci_x_W3_2T(1,:)/(dia_prtc*vg)/10^6,...
    ci_x_W3_2T(2,:)/(dia_prtc*vg)/10^6,'->b','linewidth',lin_wid);
errorbar(freq*dia_prtc/vg+0.5*delta,psd_bst_z_W3_2T/(dia_prtc*vg)/10^6,ci_z_W3_2T(1,:)/(dia_prtc*vg)/10^6,...
    ci_z_W3_2T(2,:)/(dia_prtc*vg)/10^6,':>b','linewidth',lin_wid2);
errorbar(freq*dia_prtc/vg+1.5*delta,psd_bst_x_still_2T/(dia_prtc*vg)/10^6,ci_x_still_2T(1,:)/(dia_prtc*vg)/10^6,...
    ci_x_still_2T(2,:)/(dia_prtc*vg)/10^6,'-<k','linewidth',lin_wid);
errorbar(freq*dia_prtc/vg+1.5*delta,psd_bst_z_still_2T/(dia_prtc*vg)/10^6,ci_z_still_2T(1,:)/(dia_prtc*vg)/10^6,...
    ci_z_still_2T(2,:)/(dia_prtc*vg)/10^6,':<k','linewidth',lin_wid2);

if mde==1
    ylim([0 .45])
    yticks(0:0.1:0.4);
else
    ylim([0 1.5])
end
xlabel('$fd_p/v^{\prime}_g$','interpreter','latex');
ylabel('$P(f) /v^{\prime}_g d_p$','interpreter','latex');
xlim([0 0.15])

plot(ones(2,1)*f_wave(1)*dia_prtc/vg,[0 0.5],'--r','linewidth',lin_wid)
plot(ones(2,1)*f_wave(2)*dia_prtc/vg,[0 0.5],'--g','linewidth',lin_wid)
plot(ones(2,1)*f_wave(3)*dia_prtc/vg,[0 0.5],'--b','linewidth',lin_wid)

xticks(freq*dia_prtc/vg)

grid on
xtickformat('%.3f')
legend('$(v^{\prime}-u^{\prime})^*_x$ (W1)','$(v^{\prime}-u^{\prime})^*_z$ (W1)','$(v^{\prime}-u^{\prime})^*_x$ (W2)',...
    '$(v^{\prime}-u^{\prime})^*_z$ (W2)', '$(v^{\prime}-u^{\prime})^*_x$ (W3)','$(v^{\prime}-u^{\prime})^*_z$ (W3)',...
    '$v^{\prime*}_x$ (Q)','$v^{\prime*}_z$ (Q)','interpreter','latex','numcolumns',1)

set(gcf,'position',[100 100 1000 400])

set(h(1),'position',[0.1 0.15 0.375 0.8])
set(h(2),'position',[0.575 0.15 0.375 0.8])
% set(h(3),'position',[0.1 0.075 0.375 0.375])
% set(h(4),'position',[0.575 0.075 0.375 0.375])

set(findall(h,'-property','FontSize'),'FontSize',12)
text(h(1),-0.025,0.45,'(a)','interpreter','latex','FontSize',15)
text(h(2),-0.025,0.45,'(b)','interpreter','latex','FontSize',15)
%% Representative images for N2 particle settling in quiescent fluid and waves 
% figure 12 and 13

clc
close all
load N2_spectra_25_raw_data % raw data of N2 particle in quiescent fluid and all wave cases.

% Each particle dataset (e.g., nylon_332, nylon_18, etc.) contains 10 experiments (cells).
% Column structure:
% 1st: time [s]
% 2nd-3rd: position x, z [mm]
% 4th-5th: particle velocity v_x, v_z [mm/s]
% 6th-7th: undisturbed flow velocity u_x, u_z at particle position [mm/s]

%
lin_wid1 = 0.5;
px_mm = 7.14;
z_surf = -100;
fs=20;
N2_vs = 94.6 %mm/s
freq_min = 0.606;
freq_max = freq_min*2;
T_hlf_min = 1/freq_max/2;
T_hlf_max = 1/freq_min/2;
h_dep=315;

Q_test = [14 7 9 20];
tt_auto=(0:size(W1_nylon_18_x,1)-1)/fs;

%% N2 still water spectra
clc
close all

for kk=1:length(Q_test)
    clear xx
    xx = still_nylon_18_x(:,Q_test(kk));
    xx=xx-mean(xx); 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;
    psd_test_x_still_nomean(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end

for kk=1:length(Q_test)
    clear xx
    xx = still_nylon_18_x(:,Q_test(kk));
    xx=xx; 
    num_sam = length(xx);
    test_x = xx(1:num_sam);
    xdft = fft(test_x);
    xdft = xdft(1:num_sam/2+1);
    psdx = (1/(fs*num_sam)) * abs(xdft).^2;             
    psd_test_x_still_yesmean(:,kk) = 2*psdx(1:end);
    freq = 0:fs/num_sam:fs/2;
end


%% figure 12: N2 particle settling in quiescent fluid
clc
close all

h(1) = subplot(151);
hold on
for kk=1:length(Q_test)
    test_kk=Q_test(kk);
    if kk<3
        plot((filt_nylon_18_still{test_kk}(:,2)-filt_nylon_18_still{test_kk}(1,2))/(px_mm*h_dep),...
        (filt_nylon_18_still{test_kk}(:,3)-z_surf)/(px_mm*h_dep),'-o');
    else
        plot((filt_nylon_18_still{test_kk}(:,2)-filt_nylon_18_still{test_kk}(1,2))/(px_mm*h_dep),...
        (filt_nylon_18_still{test_kk}(:,3)-z_surf)/(px_mm*h_dep),'-*');
    end
end
box on
ylim([-0.8 -0.3])
xlim([-0.1 0.1])
xlabel('$x/h$','interpreter','latex')
ylabel('$z/h$','interpreter','latex')
h(2) = subplot(152);
hold on
for kk=1:length(Q_test)
    if kk<3
        plot(tt_auto,still_nylon_18_x(:,Q_test(kk))/N2_vs,'-o');
    else
        plot(tt_auto,still_nylon_18_x(:,Q_test(kk))/N2_vs,'-*');
    end
end
box on
xlim([0 1.6])
grid on
ylim([-0.3 0.3])
yticks([-0.3:0.1:0.3])
xlabel('$t^{\prime}$ (s)','interpreter','latex')
ylabel('$v^{\prime}_x/v^{\prime}_s$','interpreter','latex')

h(3)=subplot(153)
hold on
for kk=1:length(Q_test)
    if kk<3
        plot(tt_auto,auto_x_still(Q_test(kk),:),'-o');
    else
        plot(tt_auto,auto_x_still(Q_test(kk),:),'-*');
    end
end
grid on
ylim([-1 1])
xlim([0 1.6])
xlabel('$\textrm{lag}$ (s)','interpreter','latex')
ylabel('$R_{xx} (\textrm{lag})$','interpreter','latex')
box on


h(4) = subplot(154)

hold on
for kk=1:length(Q_test)
    if kk<3
        plot(freq*dia_prtc/vg,psd_test_x_still_nomean(:,kk)/(dia_prtc*vg)/10^6,'-o');
    else
        plot(freq*dia_prtc/vg,psd_test_x_still_nomean(:,kk)/(dia_prtc*vg)/10^6,'-*');
    end
end
box on

xlim([0 0.15])
xticks(freq*dia_prtc/vg)
grid on
xlabel('$fd_p/v^{\prime}_g$','interpreter','latex');
ylabel('$P(f) /v^{\prime}_g d_p$','interpreter','latex');
box on
xtickformat('%.3f')


set(gcf,'position',[100 100 1400 600])

set(h(1),'position',[0.07 0.125 0.25 0.825])
set(h(2),'position',[0.38 0.6 0.25 0.35])
set(h(3),'position',[0.38 0.125 0.25 0.35])
set(h(4),'position',[0.7 0.6 0.25 0.35])

set(findall(h,'-property','FontSize'),'FontSize',14)
text(h(1),-0.145,-0.3,'(a)','interpreter','latex','FontSize',18)
text(h(2),-0.33,0.3,'(b)','interpreter','latex','FontSize',18)
text(h(3),-0.33,1,'(c)','interpreter','latex','FontSize',18)
text(h(4),-0.032,0.6,'(d)','interpreter','latex','FontSize',18)

legend(h(1),'Oscillating oblique','Oscillating oblique','Chaotic','Chaotic',...
    'interpreter','latex','position',[0.72 0.22 0.2 0.14],'Fontsize',16)

% (optional) save image
% saveas(figure(1),'representative_N2_Q.fig')
% saveas(figure(1),'representative_N2_Q.jpg')
% saveas(figure(1),'representative_N2_Q.eps','epsc')

%% figure 13: N2 particle settling in waves
clc
close all
% first 2 cases: oblique oscillating
% last 2 cases: chaotic
W1_test = [15 8 1 17];
W2_test = [20 17 1 2];
W3_test = [1 25 13 12];

clc
close all

h(1) = subplot(231)
hold on
for kk=1:length(W1_test)
    if kk<3
        plot(tt_auto,W1_nylon_18_x(:,W1_test(kk))/N2_vs,'-o')
    else
        plot(tt_auto,W1_nylon_18_x(:,W1_test(kk))/N2_vs,'-*')
    end
end
box on
% ylim([-0.03 0.03]);
xlim([0 1.6])
grid on
ylim([-0.3 0.3])
yticks([-0.3:0.1:0.3])
xlabel('$t^{\prime}$ (s)','interpreter','latex')
ylabel('$(v^{\prime}_x-u^{\prime}_x)/v^{\prime}_s$','interpreter','latex')

h(2) = subplot(232)
hold on
for kk=1:length(W2_test)
    if kk<3
        plot(tt_auto,W2_nylon_18_x(:,W2_test(kk))/N2_vs,'-o')
    else
        plot(tt_auto,W2_nylon_18_x(:,W2_test(kk))/N2_vs,'-*')
    end
end
box on
% ylim([-0.03 0.03]);
xlim([0 1.6])
grid on
ylim([-0.3 0.3])
yticks([-0.3:0.1:0.3])
xlabel('$t^{\prime}$ (s)','interpreter','latex')
ylabel('$(v^{\prime}_x-u^{\prime}_x)/v^{\prime}_s$','interpreter','latex')

h(3) = subplot(233)
hold on
for kk=1:length(W3_test)
    if kk<3
        plot(tt_auto,W3_nylon_18_x(:,W3_test(kk))/N2_vs,'-o')
    else
        plot(tt_auto,W3_nylon_18_x(:,W3_test(kk))/N2_vs,'-*')
    end
end
box on
% ylim([-0.03 0.03]);
xlim([0 1.6])
grid on
ylim([-0.3 0.3])
yticks([-0.3:0.1:0.3])
xlabel('$t^{\prime}$ (s)','interpreter','latex')
ylabel('$(v^{\prime}_x-u^{\prime}_x)/v^{\prime}_s$','interpreter','latex')

h(4) = subplot(234)
hold on

for kk=1:length(W1_test)
    if kk<3
        plot(tt_auto,auto_x_W1(W1_test(kk),:),'-o')
    else
        plot(tt_auto,auto_x_W1(W1_test(kk),:),'-*')
    end
end
box on
xlabel('$\textrm{lag}$ (s)','interpreter','latex')
ylabel('$R_{xx} (\textrm{lag})$','interpreter','latex')
ylim([-1 1])
grid on
xlim([0 1.6])

h(5) = subplot(235)
hold on

for kk=1:length(W2_test)
    if kk<3
        plot(tt_auto,auto_x_W2(W2_test(kk),:),'-o')
    else
        plot(tt_auto,auto_x_W2(W2_test(kk),:),'-*')
    end
end
box on
xlabel('$\textrm{lag}$ (s)','interpreter','latex')
ylabel('$R_{xx} (\textrm{lag})$','interpreter','latex')
ylim([-1 1])
grid on
xlim([0 1.6])

h(6) = subplot(236)
hold on

for kk=1:length(W3_test)
    if kk<3
        plot(tt_auto,auto_x_W3(W3_test(kk),:),'-o')
    else
        plot(tt_auto,auto_x_W3(W3_test(kk),:),'-*')
    end
end
box on
xlabel('$\textrm{lag}$ (s)','interpreter','latex')
ylabel('$R_{xx} (\textrm{lag})$','interpreter','latex')
ylim([-1 1])
grid on
xlim([0 1.6])
set(gcf,'position',[100 100 1400 600])

set(findall(h,'-property','FontSize'),'FontSize',12)
text(h(1),-0.4,0.3,'(a)','interpreter','latex','FontSize',15)
text(h(2),-0.4,0.3,'(b)','interpreter','latex','FontSize',15)
text(h(3),-0.4,0.3,'(c)','interpreter','latex','FontSize',15)

text(h(4),-0.4,1,'(d)','interpreter','latex','FontSize',15)
text(h(5),-0.4,1,'(e)','interpreter','latex','FontSize',15)
text(h(6),-0.4,1,'(f)','interpreter','latex','FontSize',15)
% (optional) save image
% saveas(figure(1),'representative_N2_slip.fig')
% saveas(figure(1),'representative_N2_slip.jpg')
% saveas(figure(1),'representative_N2_slip.eps','epsc')