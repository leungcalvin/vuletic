function [mw,bw,mwerr,bwerr,epsw,epswerr,rcs,pval] = kingtest(x,y,dx,dy,fx)
% Returns a linear fit to y = m x + b + eps * f(x) with parameter uncertainties from 2d data with 2d
% error bars. 
% For best accuracy, use the dataset with smaller uncertainties as the x data.
X = [x'.^1,x'.^0,fx'];
Y = y';

B = X\Y;
m = B(1); % good ansatz
b = B(2); % good ansatz
eps=B(3); % assume we're close to King-linear

%now, let's weight by the inverse variance.
W = diag((m^2.* dx.^2 + dy.^2).^-1);
BB=(X'*W*X)\(X'*W*Y); %weighted LSQ
BBerr = (X'*W*X).^-1;
mw = BB(1);
bw = BB(2);
epsw=BB(3);
mwerr = sqrt(BBerr(1,1));
bwerr = sqrt(BBerr(2,2));
epswerr= sqrt(BBerr(3,3));
Yhat = X*BB;
rcs = sum((Yhat - Y).^2 .* diag(W))/(length(diag(W)) - 2);
pval = chi2cdf(sum((Yhat - Y).^2 .* diag(W)),length(diag(W)) - 2);
figure;hold on;
errorbar(x,y,dy,dy,dx,dx);
plot(x,Yhat);
figure;
bar(x,(Yhat - Y).*W)

end

