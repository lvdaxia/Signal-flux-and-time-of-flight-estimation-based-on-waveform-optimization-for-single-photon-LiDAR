clc;clear;close all;

t=0:0.025:(681+400)*0.025;
load('E:\12.ʵ������(��ɾ)\ʵ��3(�߾��ȵ���)\20240115\20240122-1.mat');
[num_t,num]=his_fig_plot(data_mat(:,2));
% 
% figure(1);
% plot(num_t,num);hold on;
[value,Index_num]=max(num);
zhixin1=sum(num(Index_num-7/0.025:Index_num+10/0.025).*num_t(Index_num-7/0.025:Index_num+10/0.025))/sum(num(Index_num-7/0.025:Index_num+10/0.025));

[st,st1,st2,st3,st4]=multi_guass_5(t);
% [ Vlaue_max,Index_max]=max(st);
zhixin2=sum(t.*st)/sum(st);
tran=zhixin1-zhixin2;
% figure(1);
% figure_polish('�ɼ�����', 'northeast', num_t-tran, num, 'k-', 3);
% hold on;
% figure(1);
% xlabel('');ylabel('');

% ���������״
figure_polish('�����״', 'northeast', t, st, 'r', 2);
hold on;
% ��������Ӳ�
figure_polish('����Ӳ�', 'northeast', t, st1, 'b--', 2);
figure_polish('����Ӳ�', 'northeast', t, st2, 'b--', 2);
figure_polish('����Ӳ�', 'northeast', t, st3, 'b--', 2);
% ���ư�������ͼ��������
figure_polish({'�ɼ�����','�����״', '����Ӳ�'}, 'northeast', t, st4, 'b--', 2);

xlim([5,27]);
xlabel('ʱ��/ns');ylabel('��ֵ');



function [num_t,num]=his_fig_plot(scatter_y)
    y_max=max(scatter_y+1);
    num=zeros(y_max,1);
    for i=1:length(scatter_y)
        num(scatter_y(i))=num(scatter_y(i))+1;
    end
    num_t=0.025*(1:y_max)';
%     figure(1);
%     subplot(2,1,2);
%     bar(num_t,num);title('��״ͼģʽ');xlabel('�������ʱ��/ns');ylabel('count');set(0,'defaultfigurecolor','w');grid on;
end
function [s_t,s_t1,s_t2,s_t3,s_t4]=multi_guass_5(t)
% ����̽���� 
    a1 =    0.003965;    a2 =   0.0023;   a3 =    0.002244;   a4= 0.0003295;
    b1 =     219.1-205; b2 =   220.2-205 ; b3 =  218.2-205; b4=221-205;
    c1 =   1.466; c2 =  2.074; c3 = 1.036 ; c4= 3.467 ; 

    y1=a1*exp(-((t-b1)/c1).^2);
    y2=a2*exp(-((t-b2)/c2).^2);
    y3=a3*exp(-((t-b3)/c3).^2);
    y4=a4*exp(-((t-b4)/c4).^2);
    y=y1+y2+y3+y4;
    zhixin1=sum(y.*t)/sum(y);% ��������ģ������λ��  ����ģ�ͼ����

    Ns1=a1*c1*pi^0.5;Ts1=b1;theta_s1=c1/1.414213562373095;
    Ns2=a2*c2*pi^0.5;Ts2=b2;theta_s2=c2/1.414213562373095;
    Ns3=a3*c3*pi^0.5;Ts3=b3;theta_s3=c3/1.414213562373095;
   Ns4=a4*c4*pi^0.5;Ts4=b4;theta_s4=c4/1.414213562373095;

    s_t1=Ns1/((2*pi)^0.5*theta_s1)*exp(-(t-Ts1).^2/(2*theta_s1^2));
    s_t2=Ns2/((2*pi)^0.5*theta_s2)*exp(-(t-Ts2).^2/(2*theta_s2^2));
    s_t3=Ns3/((2*pi)^0.5*theta_s3)*exp(-(t-Ts3).^2/(2*theta_s3^2));
    s_t4=Ns4/((2*pi)^0.5*theta_s4)*exp(-(t-Ts4).^2/(2*theta_s4^2));
    s_t=s_t1+s_t2+s_t3+s_t4;
end

function figure_polish(legend_content,legend_position,x,y,color,linewidth)
    plot(x,y,color,'Linewidth',linewidth);hold on;
    grid on;
    set(gcf, 'Position', [300, 300, 400, 300]); % ���磬���½�����Ϊ (300, 300)�����Ϊ 800���߶�Ϊ 600
    h_legend = legend(legend_content, 'Location', legend_position);
    % ����ͼ���ı��������СΪ 20
    set(h_legend, 'FontSize', 10);
    ax=gca;
    % ��������� 'Box' ��������Ϊ 'on'������ʾ���±߿�
    ax.Box = 'on';
    ax.LineWidth=1.5;% �ڱ߿��ϸ
    ax.YAxis.LineWidth=1.5;% ��߿��ϸ
    ax.XAxis.LineWidth=1.5;% ��߿��ϸ
    ax.GridAlpha=0.1;
    set(gcf, 'color', 'white');% ��ͼ�εı�����ɫ����Ϊ��ɫ
    % ax.XMinorGrid='on';
    % ax.YMinorGrid='on';
end