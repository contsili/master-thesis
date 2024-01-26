%% ctf275 vs ctf 151

%% Method 1
headmodel.r = 0.10;
headmodel.o = [0 0 0.04];
headmodel.unit = 'm';

load C:\Users\user\Documents\MATLAB\matlab_toolboxes\fieldtrip\fieldtrip\template\gradiometer\ctf151.mat
load C:\Users\user\Documents\MATLAB\matlab_toolboxes\fieldtrip\fieldtrip\template\gradiometer\ctf275.mat

dippos = [0 0 0.08];
dipmom = [1 0 0]';

lf151 = ft_compute_leadfield(dippos, ctf151, headmodel) * dipmom; % Alternative: ft_prepare_leadfield which is a high-level function
lf275 = ft_compute_leadfield(dippos, ctf275, headmodel) * dipmom;

sel151 = startsWith(ctf151.label, 'M'); % disregard the reference sensors
sel275 = startsWith(ctf275.label, 'M');

lf151 = lf151(sel151,:);
lf275 = lf275(sel275, :);



%%

figure
ft_plot_topo3d(ctf151.chanpos(sel151,:), lf151)
figure
ft_plot_topo3d(ctf275.chanpos(sel275,:), lf275)

% Conclusion: ctf275 has more sensor density and leads to higher measured
% values (hence darker red and blue in the topo3d)



%% Method 2a (eq. 4): ||L_OPM||/||L_SQUID|| since N_ctf275=N_ctf151

% I know from eq. 4 that std_J_ctf275 / std_J_ctf151 = ||pinv(L_ctf275)|| * N_ctf275 / ||pinv(L_ctf151)|| * N_ctf151. 
% In methods 2a and 2b I check if this relation holds


norm(lf275) / norm(lf151) 
norm(pinv(lf151)) / norm(pinv(lf275))

% the two above give the same result. This happens as norm(lf275) = 1/norm(pinv(lf275))
% norm(lf275)
% 
% ans =
% 
%    2.2512e-05
% 
% 1/norm(pinv(lf275))
% 
% ans =
% 
%    2.2512e-05


% Conclusion: ctf275 better than ctf151 as expected. 

% Now we check if ||L|| ~ sqrt(number_sensors):
sqrt(275 / 151) % close but not the same as norm(lf275) / norm(lf151)



%% Method 2b - σ_J_OPM/σ_J_SQUID

noiselevel = 1e-6;

ntime = 100000;
s = ones(1,ntime); % this is the time course of the source and here is a random selection. In NI2 course we used ni2_activation() to create the time course.

dat151 = lf151 * s + noiselevel * randn(151,ntime);
dat275 = lf275 * s + noiselevel * randn(275,ntime);

% Way1 - faster (no loop)
est151 = pinv(lf151) * dat151; % Y=Ls+N => s_estimated = pinv(L)*Y
est275 = pinv(lf275) * dat275;

% Way2
% for i=1:ntime
%   est151(i,:) = pinv(lf151) * dat151(i,:);
%   est275(i,:) = pinv(lf275) * dat275(i,:);
% end

figure
plot(s);
hold on
plot(est151);
plot(est275);

std(est151) / std(est275) % this should get closer to norm(lf275) / norm(lf151) and norm(pinv(lf151)) / norm(pinv(lf275)) as ntime increases. It never gets equal to norm(lf275) / norm(lf151) because of the randn()

% Conclusion: eq. 4 holds! And also I saw that ||pinv(L)|| = 1 / ||L||, so
% I can work with ||L|| and not ||pinv(L)|| in eq. 4

%% Method 3 - by KT
% We checked that  ||L|| is analogous sqrt(number_sensors) but what is the specific
% number that makes it an equation?

% Nonlinear least square fitting

% ctf151
% Find the index of the first non-zero element in leadfield1 
k = find(lf151, 1, 'first');

% Create the vector x starting from the index k
x = (k:length(lf151))';

y = zeros(length(x),1);
for l = 1:length(x)
    y(l) =  norm(lf151(1:x(l))); % Note: in other simulations I used norm(pinv())
end
y=y(:);

% Use the fit function with the updated fitting model
fo = fitoptions('Method', 'NonlinearLeastSquares');
ft = fittype(@(a, n, x) a * x.^n, 'options', fo);

f1=fit(x,y,ft);

% Access the coefficients
coefficients = coeffvalues(f1);

% Extract the value of n
n151= coefficients(2);
a151= coefficients(1);


% ctf275
% Find the index of the first non-zero element in leadfield1 
k = find(lf275, 1, 'first');

% Create the vector x starting from the index k
x = (k:length(lf275))';

y = zeros(length(x),1);
for l = 1:length(x)
    y(l) =  norm(lf275(1:x(l)));
end
y=y(:);

% Use the fit function with the updated fitting model
fo = fitoptions('Method', 'NonlinearLeastSquares');
ft = fittype(@(a, n, x) a * x.^n, 'options', fo);

f1=fit(x,y,ft);

% Access the coefficients
coefficients = coeffvalues(f1);

% Extract the value of n
n275 = coefficients(2);
a275= coefficients(1);


%%
275^(n275) / 151^(n151)  % 1.2618
norm(lf275) / norm(lf151) % 1.4083

% Conclusion: they are not the same because they have a different 'a' (i.e, lf151 = a151/x^n151 vs f275 = a275/x^n275 )

(a275*275^(n275)) / (a151*151^(n151)) % 1.3745 (better than 275^(n275) / 151^(n151) but still not that close to norm(lf275) / norm(lf151))

% Conclusion: it is quite random, what happens I think it depends on the
% starting points??

%%
y1=a275*x.^(n275);
plot(x, y, 'o', 'DisplayName', 'Data');
hold on;
plot(x, y1, 'r-', 'DisplayName', 'Fit');


