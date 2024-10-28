tic % Start timing

%-------------------------------------------Constant

k = 1.3806505e-23;  % Boltzmann constant, unit: J/K (Joules per Kelvin)
e = 1.6e-19;        % Elementary charge, unit: C (Coulombs)
h = 6.6260755e-34;  % Planck constant, unit: J·s (Joules seconds)
m = 9.1e-31;        % electron mass，unit：kg
yueh = 1.05457266e-34;   % 单位J·s

% -----------------------------------------Data Input

me = 1.0098 * m;    % Effective mass of electron
mh = 0.8457 * m;    % Effective mass of hole

data = [300	-0.128655023	-129.5263557	1.679989926	6.244282322
325	-0.128899001	-138.3459503	1.605924109	6.893755042
350	-0.128884387	-146.3894496	1.541751253	7.589894847
375	-0.128611182	-153.6568535	1.487471356	8.332701737
400	-0.128079385	-160.148162	    1.44308442	9.122175714
425	-0.127288997	-165.8633752	1.408590444	9.958316775
450	-0.126240018	-170.802493	    1.383989429	10.84112492
475	-0.124932447	-174.9655155	1.369281373	11.77060016
500	-0.123366285	-178.3524425	1.364466278	12.74674247
525	-0.121541531	-180.9632743	1.369544144	13.76955188
550	-0.119458186	-182.7980106	1.384514969	14.83902837
575	-0.11711625	    -183.8566516	1.409378755	15.95517194
600	-0.114515722	-184.1391972	1.444135501	17.1179826
625	-0.111656603	-183.6456475	1.488785208	18.32746035
650	-0.108538892	-182.3760024	1.543327874	19.58360518
675	-0.10516259	    -180.3302619	1.607763501	20.8864171
700	-0.101527697	-177.5084261	1.682092089	22.2358961
725	-0.097634212	-173.9104949	1.766313636	23.63204219
];   % Temperature (T); Hall coefficient (m^3/C); Seebeck coefficient (μV K^-1); Thermal conductivity (W m^-1 K^-1); Electrical resistivity (μohm m)"

% -----------------------------------------Data reading

T = data(:,1);  % Temperature
S = data(:,3);  % Seebeck coefficient
r = data(:,5);  % Electrical resistivity
sigma = 1e-5.*1e6./r;   % Electrical conductivity, unit: 1e5 S m^-1
RH_hall = data(:,2); %  Hall coefficient
n = abs(1./RH_hall./e)./1e19;  % Carrier concentration, unit: 10^19 cm^-3

% -----------------------------------------Set calculation accuracy

p = length(T);  
ne_min= n./2;   % Minimum electron concentration
ne_max = n;  % Maximum electron concentration
step = 100;  % Appropriately increasing the step size can improve accuracy，recommended to be greater than 10


% ----------------------------------------Generate parameters
for i=1:p
    j=1:step;
    ne(i,j) = linspace(ne_min(i), ne_max(i), step); 
end

% ----------------------------------------Generate initial parameters of the same matrix size

S_expand = repmat(S,1,step);  
sigma_expand = repmat(sigma,1,step);
n_expand = repmat(n,1,step);
T_expand = repmat(T,1,step);

% ----------------------------------------Hole concentration

nh = n_expand-ne;       

% ----------------------------------------Solve yitae

f =@(yitae)(4*pi*(2*k.*me.*T_expand).^(3/2)/h.^3.*FDI(1/2,yitae).*1e-25)-ne;
yitae_guess = zeros(size(T_expand)); 
yitae = fsolve(f, yitae_guess);
% disp('找到的零点 yitae = ');
% disp(yitae);

% ----------------------------------------Solve yitah

f =@(yitah)(4*pi*(2*k.*mh.*T_expand).^(3/2)/h.^3.*FDI(1/2,yitah).*1e-25)-nh;
yitah_guess = zeros(size(T_expand));
yitah = fsolve(f, yitah_guess);
% disp('找到的零点 yitah = ');
% disp(yitah);

% ------------------------------------------Fermi-Dirac integral under the phonon scattering mechanism

FE0 = integral(@(x)x.^0./(1+exp(x-yitae)),0,100,'ArrayValued',true);
FE1 = 2.*integral(@(x)x.^1./(1+exp(x-yitae)),0,100,'ArrayValued',true);
FE2 = 3.*integral(@(x)x.^2./(1+exp(x-yitae)),0,100,'ArrayValued',true);
FH0 = integral(@(x)x.^0./(1+exp(x-yitah)),0,100,'ArrayValued',true);
FH1 = 2.*integral(@(x)x.^1./(1+exp(x-yitah)),0,100,'ArrayValued',true);
FH2 = 3.*integral(@(x)x.^2./(1+exp(x-yitah)),0,100,'ArrayValued',true);

% ------------------------------------------Calculate Se and Sh

Se = -(k/e)*(FE1./FE0-yitae)*1e6;   % Unit: μV/K 
Sh = (k/e)*(FH1./FH0-yitah)*1e6;    % Unit: μV/K 

% ------------------------------------------Calculate sigmae and sigmah

f = @(sigmae)(((Se.*sigmae+Sh.*(sigma_expand-sigmae))./sigma_expand))-S_expand;
sigmae_guess = zeros(size(T_expand)); 
options = optimoptions('fsolve', 'TolFun', 1e-1, 'TolX', 1e-1);
sigmae = fsolve(f, sigmae_guess, options);
% disp('sigmae = ');
% disp(sigmae);
sigmah = sigma_expand-sigmae;  
% disp('sigmah = ');
% disp(sigmah);

% ------------------------------------------Mobility calculated from the electrical conductivity formula

miue_sigmae = (sigmae*1e5)./(ne.*1e25)./e.*1e4;   % Unit: 10^4 cm² V^-1 s^-1
miuh_sigmah = (sigmah*1e5)./(nh.*1e25)./e.*1e4;   % Unit: 10^4 cm² V^-1 s^-1

% ------------------------------------------Calculate RH using the two-band model

RH = (nh.*miuh_sigmah.^2-ne.*miue_sigmae.^2)./(ne.*miue_sigmae+nh.*miuh_sigmah).^2./e./1e19; 
RH_hall_expand = repmat(RH_hall,1,step);  % Experimental measurement of the Hall coefficient

% ------------------------------------------Calculate relative error

err_relative = RH_hall_expand-RH;
err = abs(RH_hall_expand - RH)./abs(RH_hall)*100;

% ------------------------------------------Calculate the value of the entire matrix A

A =me^1.5.*miue_sigmae./mh^1.5./miuh_sigmah;

% ------------------------------------------Initialize the position parameter of the minimum error

min_error_positions = zeros(size(err, 1), 1);
min_errors = zeros(size(err, 1), 1);

% ------------------------------------------Output the position parameter of the minimum error

for i = 1:size(err, 1)
    % 获取第i行的最小值及其位置
    [min_value, min_index] = min(err(i, :));
    min_error_positions(i) = min_index;
    min_errors(i) = min_value;
end

% ------------------------------------------Output other parameters corresponding to that position

for i = 1:size(err, 1)
   disp(['Row', num2str(i), 'minimum value position: ', num2str(min_error_positions(i))]);
   disp(['Row', num2str(i), 'minimum error: ', num2str(min_errors(i))]); 
   best_sigmah(i) = sigmah (i,min_error_positions(i));
   best_sigmae(i) = sigmae (i,min_error_positions(i));
   best_ne(i) = ne(i,min_error_positions(i));
   best_nh(i) = nh(i,min_error_positions(i));
   best_Se(i) = Se(i,min_error_positions(i));
   best_Sh(i) = Sh(i,min_error_positions(i));
   best_miue(i) = miue_sigmae(i,min_error_positions(i));
   best_miuh(i) = miuh_sigmah(i,min_error_positions(i));
   best_RH(i) = RH(i,min_error_positions(i));
   best_FE0(i) = FE0(i,min_error_positions(i));
   best_FE1(i) = FE1(i,min_error_positions(i));
   best_FE2(i) = FE2(i,min_error_positions(i));
   best_FH0(i) = FE0(i,min_error_positions(i));
   best_FH1(i) = FE1(i,min_error_positions(i));
   best_FH2(i) = FH2(i,min_error_positions(i));
   best_A(i)= A(i,min_error_positions(i));
end

% ------------------------------------------Calculate other parameters

 T =T';
 Kbip = (best_sigmae.*best_sigmah./(best_sigmae+best_sigmah)).*(best_Se-best_Sh).^2.*T./1e7; % Bipolar thermal conductivity
 Le = (k/e)^2.*(best_FE2./best_FE0-(best_FE1./best_FE0).^2); % Electronic contribution Lorenz number
 Lh = (k/e)^2.*(best_FH2./best_FH0-(best_FH1./best_FH0).^2); % Hole contribution Lorenz number
 Kele = (Le.*best_sigmae.*T+Lh.*best_sigmah.*T).*1e5; % Electronic thermal conductivity
 Klat = data(:,4)-Kele'-Kbip'; % Lattice thermal conductivity
 
 % ------------------------------------------Output all parameters
 
 output = [min_error_positions,best_RH',RH_hall,min_errors,T',best_ne',best_nh',best_Se',best_Sh',best_miue',best_miuh',best_sigmae',best_sigmah',Kbip',Kele',Klat,best_A'];

 toc % End timing