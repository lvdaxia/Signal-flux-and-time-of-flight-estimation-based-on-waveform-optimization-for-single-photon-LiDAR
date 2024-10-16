clc;clear;close all;

t=0:0.025:(681+400)*0.025;
load('E:\12.实验数据(勿删)\实验3(高精度调试)\20240115\20240122-1.mat');
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
% figure_polish('采集波形', 'northeast', num_t-tran, num, 'k-', 3);
% hold on;
% figure(1);
% xlabel('');ylabel('');

% 绘制拟合形状
figure_polish('拟合形状', 'northeast', t, st, 'r', 2);
hold on;
% 绘制拟合子波
figure_polish('拟合子波', 'northeast', t, st1, 'b--', 2);
figure_polish('拟合子波', 'northeast', t, st2, 'b--', 2);
figure_polish('拟合子波', 'northeast', t, st3, 'b--', 2);
% 绘制包含两个图例的曲线
figure_polish({'采集波形','拟合形状', '拟合子波'}, 'northeast', t, st4, 'b--', 2);

xlim([5,27]);
xlabel('时间/ns');ylabel('幅值');



function [num_t,num]=his_fig_plot(scatter_y)
    y_max=max(scatter_y+1);
    num=zeros(y_max,1);
    for i=1:length(scatter_y)
        num(scatter_y(i))=num(scatter_y(i))+1;
    end
    num_t=0.025*(1:y_max)';
%     figure(1);
%     subplot(2,1,2);
%     bar(num_t,num);title('柱状图模式');xlabel('脉冲飞行时间/ns');ylabel('count');set(0,'defaultfigurecolor','w');grid on;
end
function [s_t,s_t1,s_t2,s_t3,s_t4]=multi_guass_5(t)
% 线性探测器 
    a1 =    0.003965;    a2 =   0.0023;   a3 =    0.002244;   a4= 0.0003295;
    b1 =     219.1-205; b2 =   220.2-205 ; b3 =  218.2-205; b4=221-205;
    c1 =   1.466; c2 =  2.074; c3 = 1.036 ; c4= 3.467 ; 

    y1=a1*exp(-((t-b1)/c1).^2);
    y2=a2*exp(-((t-b2)/c2).^2);
    y3=a3*exp(-((t-b3)/c3).^2);
    y4=a4*exp(-((t-b4)/c4).^2);
    y=y1+y2+y3+y4;
    zhixin1=sum(y.*t)/sum(y);% 激光脉冲模型质心位置  理论模型计算的

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
    set(gcf, 'Position', [300, 300, 400, 300]); % 例如，左下角坐标为 (300, 300)，宽度为 800，高度为 600
    h_legend = legend(legend_content, 'Location', legend_position);
    % 设置图例文本的字体大小为 20
    set(h_legend, 'FontSize', 10);
    ax=gca;
    % 将坐标轴的 'Box' 属性设置为 'on'，以显示上下边框
    ax.Box = 'on';
    ax.LineWidth=1.5;% 内边框粗细
    ax.YAxis.LineWidth=1.5;% 外边框粗细
    ax.XAxis.LineWidth=1.5;% 外边框粗细
    ax.GridAlpha=0.1;
    set(gcf, 'color', 'white');% 将图形的背景颜色设置为白色
    % ax.XMinorGrid='on';
    % ax.YMinorGrid='on';
end