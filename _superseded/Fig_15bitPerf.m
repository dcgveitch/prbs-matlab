figure
hold on
colours=pmkmp(6,'CubicL');
X=[2 3 4 6 12 24];

for d_i=unique(out_summary(:,1))'
    out_plot{d_i}=out_summary(find(out_summary(:,1)==d_i & out_summary(:,3)==2),[11 13 15]);
    plot(X,out_plot{d_i}(:,1),'.--','color',colours(d_i+1,:));
    plot(X,out_plot{d_i}(:,2),'.-','color',colours(d_i+1,:));
    plot(X,out_plot{d_i}(:,3),'.--','color',colours(d_i+1,:));
end
ylim([-0.5 0.5]);
export_fig('filename', fileSaveName2, '-nocrop');
close(gcf);