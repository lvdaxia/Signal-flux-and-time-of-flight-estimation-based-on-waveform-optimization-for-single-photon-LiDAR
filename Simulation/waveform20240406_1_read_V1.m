% 本函数用于查看拟合参数，方差等
% waveform20240405_4读取
% 时间分辨力25ps，累积时间10min*60s=600s，重频3000Hz
addpath('D:\matlab\bin\SPAD\muiti_guss\TIM_measurement_V2\waveform_characterization');
clc;clear;%close all;
load('E:\12.实验数据(勿删)\实验3(高精度调试)\20240406\20240405-3k-waveform-4.mat');
chongpin=3000;
Time_length=0.1*round(data_mat(end,1)/1e9/3600*10)*3600;
upper=228;
lower=210;
figure(1);
scatter(data_mat(:,1)/1e9/3600,data_mat(:,2)*0.025,'.');xlabel('时间/小时h');ylabel('飞行时间/ns');
scatter_x=data_mat(:,1);
scatter_y=data_mat(:,2);
% 
% scatter_x=scatter_x(scatter_y*0.025<228);
% scatter_y=scatter_y(scatter_y*0.025<228);

[num_t,num]=his_fig_plot(scatter_y);
figure(6);
plot(num_t,num);grid on;hold on;
title('剪切位移前的波形');xlabel('时间/ns');ylabel('计数/counts');
P=sum(num)/(600*3000);%10分钟
Ns_solve=-log(1-P);
fprintf('测量时间：%4d\n',data_mat(end,1)/1e9);
fprintf('信号光子数：%1.4d\n',Ns_solve);
time_ac=ceil(data_mat(end,1)/1e9*3000)/3000;
[Ns_recover_boost,Ns_recover,Index_upper,Index_lower]=correction(num, time_ac);
Ns_recover_boost_x=num_t(Index_lower:Index_upper);
zhixin_correct=sum(Ns_recover_boost_x.*Ns_recover_boost)/sum(Ns_recover_boost);
figure(6);
plot(Ns_recover_boost_x,Ns_recover_boost,'r');hold on;
plot(Ns_recover_boost_x,Ns_recover*time_ac*3000,'g');hold on;

% figure(6);
% Ns_recover_nor=Ns_recover_boost/sum(Ns_recover_boost);
% tra_recover_nor=Ns_recover_boost_x-(zhixin_correct-219.5095);
% plot(tra_recover_nor,Ns_recover_nor,'r');hold on;
xlabel('time of flight/ns');ylabel('counts');
figure_polish();
% [fitresult, gof] = createFit(tra_recover_nor, Ns_recover_nor);
[tw]=Tw_calculation_f(Ns_recover_boost_x, Ns_recover_boost);


function figure_polish()
    grid on;
    set(gcf, 'Position', [300, 300, 400, 300]); % 例如，左下角坐标为 (300, 300)，宽度为 800，高度为 600
    %h_legend = legend(legend_content, 'Location', legend_position);
    % 设置图例文本的字体大小为 20
%     set(h_legend, 'FontSize', 10);
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

function [num_t,num]=his_fig_plot(scatter_y)
    y_max=max(scatter_y);
    num=zeros(y_max,1);
    for i=1:length(scatter_y)
        num(scatter_y(i))=num(scatter_y(i))+1;
    end
    num_t=(1:y_max)*0.025;
end


function [Ns_recover_boost,Ns_recover,Index_upper,Index_lower]=correction(num, time_ac)
    [~,Index_max]=max(num);

    Index_upper=Index_max+15/0.025;
    Index_lower=Index_max-7/0.025;

    laser_frequency=3000;

    for i=1:Index_upper
        P(i)=num(i)/ceil(laser_frequency*time_ac);% sum_num/(重频*累积时间)  注意设置参数
    end

    for i=Index_lower:Index_upper
%         P(i)=num(i)/ceil(laser_frequency*time_ac);% sum_num/(重频*累积时间)  注意设置参数
        Ns_recover(i-Index_lower+1)=-log(1-P(i)/(1-sum(P(i-77/0.025:i-1))));
        if isnan(Ns_recover(i-Index_lower+1))  
            Ns_recover(i-Index_lower+1)=0;
        end
        if ~isreal(Ns_recover(i-Index_lower+1))
            Ns_recover(i-Index_lower+1)=inf;
        end
    end

    Ns_recover_boost=Ns_recover*sum(num(Index_lower:Index_upper))/sum(Ns_recover);
    num_t=(Index_lower:Index_upper)*0.025;
    zhixin=sum(num_t.*num(Index_lower:Index_upper)')/sum(num(Index_lower:Index_upper));
%     figure(1);
%     plot(num_t,num);hold on;
%     plot(num_t(Index_lower:Index_upper),Ns_recover);hold on;
%     plot(num_t(Index_lower:Index_upper),Ns_recover_boost);hold on;
%     xlim([208,224]);
end
function [fitresult, gof] = createFit(tra_recover_nor, Ns_recover_nor)
    %% Fit: 'untitled fit 1'.
    [xData, yData] = prepareCurveData( tra_recover_nor, Ns_recover_nor );

    % Set up fittype and options.
    ft = fittype( 'gauss4' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf 0 0 -Inf 0 -Inf -Inf 0 -Inf -Inf 0];
    opts.Robust = 'LAR';
    opts.StartPoint = [0.00728031541640431 218.912491690585 0.332727964809446 0.00647283338006824 219.587491690585 0.369891396507293 0.00617272709588556 218.362491690585 0.464985530935887 0.00476201690142912 220.312491690585 0.51290597245554];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult, xData, yData );
    legend( h, 'Ns_recover_nor vs. tra_recover_nor', 'untitled fit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel tra_recover_nor
    ylabel Ns_recover_nor
    grid on
end
function [tw]=Tw_calculation_f(t2, s2)
    t_pin=sum(t2.*s2)/sum(s2);
    tw=(sum((t2-t_pin).^2.*s2)/sum(s2))^0.5;
end

