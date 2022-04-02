%% Load Data
% 23 x 28 x 1497 (23 trials 28 channels 1479 data)
load new_sig_close_a;% loading data
data=new_sig_close_a(1,:,:);
D=squeeze(data);
labls={'fp_1','f_7','f_3','fz','fc_5','fc_1','t_7','c_3',...%left frontal
       'cp_5','cp_1','p_7','p_3','o_1',...%
       'fp_2','f_8','f_4','fC_6','fC_2','t_8','c_4','cz',...%
       'cP_6','cp_2','p_8','p_4','Pz','o_2','oz'};%
Fs=512; % sampling frequency
T=1/Fs;
[NCh,LL]=size(D);
t=(0:LL-1)*T;    % Time vector;%generates time vector
figure(1)
title('EEG signal')
for i=1:NCh
    subplot(NCh,1,i);plot(t,D(i,:))
    if i==1;title('EEG channels');end
    xlim([0 t(end)]);ylabel([labls{i}],'FontSize',10)
    ax=gca;ax.YTick=[];if i<NCh; ax.XTick=[];end
    ylh = get(gca,'ylabel');gyl = get(ylh);ylp = get(ylh, 'Position');
    set(ylh,'Rotation',0,'Position',ylp,'VerticalAlignment','middle','HorizontalAlignment','right')
end
xlabel('Time(sec)');
%% DELTA band 1-4 Hz
Dpass = 0.057501127785;  % Passband Ripple
Dstop = 0.0001;          % Stopband Attenuation
dens  = 20;              % Density Factor
[N, Fo, Ao, W] = firpmord([0,4]/(Fs/2), [1 0], [Dpass, Dstop]);
b1 = firpm(N, Fo, Ao, W, {dens});Hd1 = dfilt.dffir(b1);
x1=filter(Hd1,D)';
% Wn=[1 4]/Fs;
% [b1,a1] = butter(2,Wn,'bandpass');
% x2 = filter(b1,a1,D)'; 
figure(2);A=subplot(521);plot(t,x1(:,1:2)','LineWidth',1.5);xlim([0 t(end)]);
legend('Ch1','Ch2','Orientation','horizontal');title('Delta 1-4Hz','FontSize',10);
pos=get(A,'Position');set(gca,'xtick',[],'position',[pos(1)-0.05,pos(2),pos(3)+0.12,pos(4)+0.025]);
L = 512; noverlap =0.5*L; [ps2,f] = pwelch(x1,hamming(L),noverlap);pdelta=10*log10(ps2);
A=subplot(522);plot(f,pdelta(:,1:2),'LineWidth',1.5);xlim([0 t(end)]/2);title('PSD (dB/Hz)','FontSize',10);
pos=get(A,'Position');set(gca,'xtick',[],'position',[pos(1)+0.05,pos(2),pos(3),pos(4)+0.025]);ylim([-150 -20])
%% Tetha band 4-8Hz
Fstop1 = 3.5;             % First Stopband Frequency
Fpass1 = 4;               % First Passband Frequency
Fpass2 = 7;               % Second Passband Frequency
Fstop2 = 7.5;             % Second Stopband Frequency
Dstop1 = 0.001;           % First Stopband Attenuation
Dpass  = 0.057501127785;  % Passband Ripple
Dstop2 = 0.0001;          % Second Stopband Attenuation
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 0], [Dstop1 Dpass Dstop2]);
b2 = firpm(N, Fo, Ao, W, {dens});Hd2 = dfilt.dffir(b2);x2=filter(Hd2,D)';
figure(2);A=subplot(523);plot(t,x2(:,1:2),'LineWidth',1.5);xlim([0 t(end)]);title('Tetha 4-8Hz','FontSize',10);
pos=get(A,'Position');set(gca,'xtick',[],'position',[pos(1)-0.05,pos(2),pos(3)+0.12,pos(4)+0.025]);
[ps2,f] = pwelch(x2,hamming(L),noverlap);ptetha=10*log10(ps2);
A=subplot(524);plot(f,ptetha(:,1:2),'LineWidth',1.5);xlim([0 t(end)]/2);
pos=get(A,'Position');set(gca,'xtick',[],'position',[pos(1)+0.05,pos(2),pos(3),pos(4)+0.025]);ylim([-150 -20])
%% ALPHA band 8-13
Fstop1 = 7.5;             % First Stopband Frequency
Fpass1 = 8;               % First Passband Frequency
Fpass2 = 12;              % Second Passband Frequency
Fstop2 = 12.5;            % Second Stopband Frequency
Dstop1 = 0.001;          % First Stopband Attenuation
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 0], [Dstop1 Dpass Dstop2]);
b3  = firpm(N, Fo, Ao, W, {dens});Hd3 = dfilt.dffir(b3); x3=filter(Hd3,D)';
figure(2); A=subplot(525);plot(t,x3(:,1:2),'LineWidth',1.5);xlim([0 t(end)]);title('Alpha 8-13Hz','FontSize',10);
pos=get(A,'Position');set(gca,'xtick',[],'position',[pos(1)-0.05,pos(2),pos(3)+0.12,pos(4)+0.03]);
[ps2,f] = pwelch(x3,hamming(L),noverlap);palpha=10*log10(ps2);
A=subplot(526);plot(f,palpha(:,1:2),'LineWidth',1.5);xlim([0 t(end)]/2);
pos=get(A,'Position');set(gca,'xtick',[],'position',[pos(1)+0.05,pos(2),pos(3),pos(4)+0.025]);ylim([-150 -20])
%% BETA band 13-30
Fstop1 = 11.5;            % First Stopband Frequency
Fpass1 = 12;              % First Passband Frequency
Fpass2 = 30;              % Second Passband Frequency
Fstop2 = 30.5;            % Second Stopband Frequency
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 0], [Dstop1 Dpass Dstop2]);
b4 = firpm(N, Fo, Ao, W, {dens});Hd4 = dfilt.dffir(b4);x4=filter(Hd4,D)';
figure(2);A=subplot(527);plot(t,x4(:,1:2),'LineWidth',1.5);xlim([0 t(end)]);title('Beta 13-30Hz','FontSize',10);
pos=get(A,'Position');set(gca,'xtick',[],'position',[pos(1)-0.05,pos(2),pos(3)+0.12,pos(4)+0.03]);

[ps2,f] = pwelch(x4,hamming(L),noverlap);pbeta=10*log10(ps2);
A=subplot(528);plot(f,pbeta(:,1:2),'LineWidth',1.5);xlim([0 t(end)]/2);
pos=get(A,'Position');set(gca,'xtick',[],'position',[pos(1)+0.05,pos(2),pos(3),pos(4)+0.025]);ylim([-150 -20])
%% Gamma band 30-45
Fstop1 = 29.5;            % First Stopband Frequency
Fpass1 = 30;              % First Passband Frequency
Fpass2 = 45;              % Second Passband Frequency
Fstop2 = 45.5;            % Second Stopband Frequency
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), [0 1 0], [Dstop1 Dpass Dstop2]);
b5 = firpm(N, Fo, Ao, W, {dens});Hd5 = dfilt.dffir(b5);x5=filter(Hd5,D)';
figure(2);A=subplot(529);plot(t,x5(:,1:2),'LineWidth',1.5);xlim([0 t(end)]);title('Gamma 30-45Hz','FontSize',10);
xlabel('Time(sec)');pos=get(A,'Position');set(gca,'xtick',[],'position',[pos(1)-0.05,pos(2),pos(3)+0.12,pos(4)+0.03]);

[ps2,f] = pwelch(x5,hamming(L),noverlap);pgamma=10*log10(ps2);
A=subplot(5,2,10);plot(f,pgamma(:,1:2),'LineWidth',1.5);xlim([0 t(end)]/2);xlabel('Frequency (Hz)');
pos=get(A,'Position');set(gca,'xtick',[],'position',[pos(1)+0.05,pos(2),pos(3),pos(4)+0.025]);ylim([-150 -20])
%% Define electrodes placement for 32 Ch
Elec=xlsread('LocatElectrode','Ch32');
LCh32={'Fp_1','AF_3','F_7','F_3','FC_1','FC_5','T_7','C_3','CP_1','CP_5','P_7',...
       'P_3','Pz','PO_3','O_1','Oz','O_2','PO_4','P_4','P_8','CP_6','CP_2','C_4',...
       'T_8','FC_6','FC_2','F_4','F_8','AF_4','Fp_2','Fz','Cz'};
x = Elec(:,1);
y = Elec(:,2);
z = Elec(:,3);
t=t';
figure(4)
xx=zeros(32,1);yy=xx;zz=xx;
labels=cell(32,1);
for i=1:32
    k=strcmpi(labls,LCh32(i));
    k(k==0)=[];
    if k==1
       plot3(x(i)-1,y(i)-3,z(i)+6,'ro','MarkerSize',1,'LineWidth',5,...
            'MarkerEdgeColor','r')
       set(gca,'color',[0.97 0.9 0.5])
       text(x(i)-1,y(i)+1,z(i),LCh32(i),'Color','k','FontSize',8,'FontWeight',...
            'bold','FontName','Times New Roman')
      xx(i,:)=x(i);
      yy(i,:)=y(i);
      zz(i,:)=z(i);
      labels{i,:}=LCh32{i};
    else
        plot3(x(i)-1,y(i)-3,z(i)+6,'bo','MarkerSize',1,'MarkerEdgeColor','b',...
              'LineWidth',5)
        text(x(i)-1.5,y(i)+1,z(i),LCh32(i),'Color','k','FontSize',8,'FontWeight',...
            'bold','FontName','Times New Roman')
      xx(i,:)=inf;
      yy(i,:)=inf;
      zz(i,:)=inf;
      labels{i,:}='inf';
    end
    xlabel('X');ylabel('Y');zlabel('Z')
    title('Electrodes position for 32 channels','FontSize',12,'FontWeight',...
          'bold','FontName','Times New Roman')
    hold  on
end
legend({'Active','Non-Active'},'FontSize',10,'FontWeight','bold','Location','best')
labels(ismember(labels,'inf'))=[];
%%
yy(yy==inf)=[];xx(xx==inf)=[];zz(zz==inf)=[];
xi=linspace(min(xx),max(xx),28);
yi=linspace(min(yy),max(yy),28);
zi=linspace(min(zz),max(zz),28);
[XI,YI]=meshgrid(xi,yi);
for j=1:5
    trlen = [];
    % Calculating power of each electrode
    if j==1
        x=x1;c='Delta (1–4 Hz)';
    elseif j==2
        x=x2;c='Theta (4–8 Hz)';
    elseif j==3
        x=x3;c='Alpha (8–13 Hz)';
    elseif j==4
        x=x4;c='Beta (13–30 Hz)';
    else
        x=x5;c='Gamma (30–45 Hz)';
    end
        for i=1:28
            ftr = x(:,i);ftr = fft(ftr);pow = ftr.*conj(ftr);
            tpow = sum(pow);trlen = vertcat(trlen,tpow);%#ok
        end
        ZI = griddata(xx,yy,trlen,XI,YI,'cubic');
        figure(5);subplot(2,5,j);contourf(XI,YI,ZI,10);hold on;
        A=scatter(xx/1.2,yy/1.2,'r','filled'); set(gca,'Visible','off');a=A.Parent.Position;
        text(a(1)-90,120,c,'FontSize',10,'FontWeight','bold')
        % %%%%%%%%%%%%%%%%%%%%%%%%%%% Draw head %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        l = 0:2*pi/100:2*pi;
        basex = .18*max(yy);  
        tip = max(yy)*1.15; base = max(yy)-.004;
        EarX = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
        EarY = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];
        plot(cos(l).*max(yy),sin(l).*max(yy),'color',[0 0 0],'Linestyle','-','LineWidth',2); % plot head
        plot([.18*max(yy);0;-.18*max(yy)],[base;tip;base],'Color',[0 0 0],'LineWidth',2); % plot nose
        plot(EarX*max(xx)+40,EarY*max(yy),'color',[0 0 0],'LineWidth',2)  % plot left ear
        plot(-EarX*max(xx)-40,EarY*max(yy),'color',[0 0 0],'LineWidth',2) % plot right ear
        box off;set(gca,'XTick',[],'YTick',[]);
        hold off
        colormap(parula);colorbar
end
%%
figure(5)
Fn = Fs/2;                             % Nyquist Frequency
Fv = linspace(0, 1, fix(LL/2)+1)*Fn;   % Frequency Vector
Iv = 1:length(Fv);                     % Index Vector
M=length(Fv);
for h=1:5
    PLI=zeros(28,28);
    if h==1
        x=x1;
    elseif h==2
        x=x2;
    elseif h==3
        x=x3;
    elseif h==4
        x=x4;
    else
        x=x5;
    end
    for i=1:28
        ffti = fft(x(:,i))/LL;            % Normalised Fourier Transform
        phi  = angle(ffti(Iv));            % Spectrum Phase *180/pi = degree
        for j=1:28
            b=0;
            fftj = fft(x(:,j))/LL;
            %amp_fts = abs(fts(Iv))*2;       % Spectrum Amplitude
            phj  = angle(fftj(Iv));
            ph=phi-phj;
            for k=1:M-1
                A=sign(ph(k));
                b=b+A;
            end
            PLI(i,j) = abs(b/M);
        end
    end
    %plotting brain network with nodes of equal sizes and no colorbars
    p = 0.07;                      %proportion of weigthed links to keep for.
    W=PLI;
    r_nodepos=[xx/1.2 yy/1.2];
    aij = threshold(W, p);         %thresholding networks due to proportion p
    ijw = adj2edgeL(triu(aij));    %passing from matrix form to edge list form
    w_atribut = ijw(:, 3);
    for lk = 1 : size(ijw, 1)      %along all links
        xynodes = zeros(2, 2);     %choose the position XY of nodes of one link
        for nd = 1 : 2             %for the two nodes of a link
            xynodes(nd, :) = r_nodepos(ijw(lk,nd), :);  %hold the positions
        end
        %  line(xynodes(:, 1), xynodes(:, 2), 'LineWidth', w_atribut_new(lk), 'Color', RGBlinks(lk, :));
         subplot(2,5,h+5)
        line(xynodes(:, 1), xynodes(:, 2), 'LineWidth', ijw(lk, 3), 'Color','r'); % links de un solo color
        hold on
        set(gca,'Visible','off');
        for nd = 1 : 2              %for the two nodes of a link
            %   plot(xynodes(nd, 1), xynodes(nd, 2), 'o', 'MarkerSize', 5*pi.*n_atribut_new(ijw(lk,nd)), ...
            %   'MarkerEdgeColor','k', 'MarkerFaceColor', RGBnodes(ijw(lk,nd), :), 'LineWidth',1.5);
            plot(xynodes(nd, 1), xynodes(nd, 2), '.k', 'MarkerSize', 2, ...
            'MarkerEdgeColor','k', 'MarkerFaceColor', 'r', 'LineWidth',1.5);
            text(xynodes(nd, 1)-0.01, xynodes(nd, 2)-0.01, labels{ijw(lk,nd)}, 'fontsize', 8);
            % text(xynodes(nd, 1)+0.01, xynodes(nd, 2)+0.01, num2str(ijw(lk,nd)), 'fontsize', 15); % number of node as label
            l = 0:2*pi/100:2*pi;
            basex = .18*max(yy);
            tip = max(yy)*1.15; base = max(yy)-.004;
            EarX = [.497 .510 .518 .5299 .5419 .54 .547 .532 .510 .489];
            EarY = [.0555 .0775 .0783 .0746 .0555 -.0055 -.0932 -.1313 -.1384 -.1199];
           
            plot(cos(l).*max(yy),sin(l).*max(yy),'color',[0 0 0],'Linestyle','-','LineWidth',2); % plot head
            plot([.19*max(yy);0;-.19*max(yy)],[base;tip;base],'Color',[0 0 0],'LineWidth',2); % plot nose
            plot(EarX*max(xx)+40,EarY*max(yy),'color',[0 0 0],'LineWidth',2)  % plot left ear
            plot(-EarX*max(xx)-40,EarY*max(yy),'color',[0 0 0],'LineWidth',2) % plot right ear
%             axis square;axis off;axis tight;box off;
            set(gca,'XTick',[],'YTick',[]);
            set(gcf, 'units','normalized','outerposition',[0 0 1 1]); %EXPANDING FIGURE ON SCREEN
        end
    end
end
%% JF JO JP JC
figure();
for h=1:5
    ph=zeros(28,length(Iv));
    if h==1
        x=x1;
        c='Delta (1–4 Hz)';
    elseif h==2
        x=x2;
        c='Theta (4–8 Hz)';
    elseif h==3
        x=x3;
        c='Alpha (8–13 Hz)';
    elseif h==4
        x=x4;
        c='Beta (13–30 Hz)';
    else
        x=x5;
        c='Gamma (30–45 Hz)';
    end
    for i=1:28
        ffti = fft(x(:,i))/LL;            % Normalised Fourier Transform
        ph(i,:)  = angle(ffti(Iv));            % Spectrum Phase *180/pi = degree
    end
    Matrix=[sum(ph(1:4,:))+sum(ph(14:16,:));ph(13,:)+sum(ph(27:28,:));...
            sum(ph(9:12,:))+sum(ph(22:26,:));sum(ph(5:8,:))+sum(ph(17:21,:))];
    PLII=zeros(4,4);
    for i=1:4
        b=0;
        for j=1:4
            g=Matrix(i,:)-Matrix(j,:);
            for k=1:M-1
                A=sign(g(k));
                b=b+A;
            end
            PLII(i,j) = abs(b/M);
        end
    end
    subplot(2,5,h)
    imagesc(PLII)
%     yticks([])
    colorbar('southoutside')
    %yticks([1 2 3 4])
%     yticklabels({'JF', 'JO', 'JP', 'JC'})
%     xticks([1 2 3 4])
%     xticklabels({'JF', 'JO', 'JP', 'JC'})
    title(num2str(c),'FontSize',16,'FontWeight','bold','FontName','Times New Roman')
end
