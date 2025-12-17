clc; clear;close all

x_num = 22; y_num = 8; src_num = x_num * y_num;
deg = 110:10:290;

load red.mat;  rgb1 = rgb;
load blue.mat; rgb2 = rgb;

fs = 12;
fontname = 'Arial';

files = {
    "src_true.txt" %test model
    "inv_result/sph1_result/src_inv_f1.txt" %inversion result
    };

labels = {'(a)','(b)','(c)','(d)'};


figure('Units','normalized','Position',[0.1 0.1 0.55 0.55]);
t = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

for i = 1:2
    da = readmatrix(files{i});
    A_init = da(:,1); A_init(isnan(A_init)) = 0;
    t_init = da(:,2); t_init(A_init<max(A_init)/1000) = 0;
    d_init = deg(da(:,3));

    A1 = flipud(reshape(A_init, [x_num, y_num])');
    t1 = flipud(reshape(t_init, [x_num, y_num])');
    d1 = flipud(reshape(d_init, [x_num, y_num])');

    ax1 = nexttile((i-1)*2+1);
    imagesc(ax1, A1); colormap(ax1, rgb1); clim(ax1, [0 10]); hold(ax1,'on');

    [cols, rows] = meshgrid(1:x_num, 1:y_num);
    dx = A1 .* cosd(d1); dy = -A1 .* sind(d1);
    quiver(ax1, cols, rows, dx, dy, 0.8, 'w', 'LineWidth', 1);

    axis(ax1,'equal','tight');
    xlim(ax1,[0.5 x_num+0.5]); ylim(ax1,[0.5 y_num+0.5]);
    set(ax1,'FontSize',fs,'FontName',fontname,'TickDir','out');
    ylabel(ax1,'Depth (km)');
    xticks(ax1,[0.5,5.5,10.5,15.5,20.5]); xticklabels(ax1,{'0','25','50','75','100'});
    yticks(ax1, [2.3475,4.6569,6.9663]); yticklabels(ax1, {'10','20','30'});
    ax1.XAxisLocation = 'top';

    ax2 = nexttile((i-1)*2+2);
    imagesc(ax2, t1); colormap(ax2, rgb2); clim(ax2,[0 200]);
    axis(ax2,'equal','tight');
    xlim(ax2,[0.5 x_num+0.5]); ylim(ax2,[0.5 y_num+0.5]);
    set(ax2,'FontSize',fs,'FontName',fontname,'TickDir','out');
    ylabel(ax2,'Depth (km)');
    xticks(ax2,[0.5,5.5,10.5,15.5,20.5]); xticklabels(ax2,{'0','25','50','75','100'});
    yticks(ax2, [2.3475,4.6569,6.9663]); yticklabels(ax2, {'10','20','30'});
    ax2.XAxisLocation = 'top';


    if i==1
        xlabel(ax1,'Distance along the strike of Fault-1 (km)');
        xlabel(ax2,'Distance along the strike of Fault-1 (km)');
    end

    if i==2
        colormap(ax1, rgb1);
        clim([0 10]);
        cb1 = colorbar(ax1,'southoutside');
        cb1.Label.String = 'Strength Coefficient';
        cb1.Label.FontSize = fs;
        cb1.FontSize = fs;
        cb1.TickDirection = 'out';

        colormap(ax2, rgb2);
        clim([0 200]);
        cb2 = colorbar(ax2,'southoutside');
        cb2.Label.String = 'Excitation Time (s)';
        cb2.Label.FontSize = fs;
        cb2.FontSize = fs;
        cb2.TickDirection = 'out';
    end

    text(ax1,-0.13,1.15,labels{(i-1)*2+1},'Units','normalized','FontName',fontname,'FontSize',fs+1);
    text(ax2,-0.13,1.15,labels{(i-1)*2+2},'Units','normalized','FontName',fontname,'FontSize',fs+1);
end



