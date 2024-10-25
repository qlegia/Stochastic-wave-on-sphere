Lmax = 600;
Nhat = 100;  % number of samples

myd= [0:0.01:pi];
UXs = zeros(Nhat,1);
UYs = zeros(Nhat, length(myd));
for inst=0:Nhat-1
    fname = sprintf('Values_Lmax%d_instance%d.mat',Lmax,inst);
    eval(['load ' fname]);
    UXs(inst+1) = Ux;
    UYs(inst+1,:) = Uy;
    mdd = dd;
end	
% compute the variance at d_k over Nhat samples
for k=1:length(dd)
    vv(k) = var(UXs(:) - UYs(:,k)); 	
end
figure('Renderer', 'painters', 'Position', [10 10 900 600]);
plot(dd,vv,dd,2.6*dd.^0.25);
xlim([0 max(dd)])
