f = 2*pi*linspace(1e5, 0.8e6, 10000);

tau = 1/1.8e6;
tau2 =1/4e4;

mu_init = 4300 ./ (1 + 1j* tau*f); 

mu_p = real(mu_init);
mu_pp = (f*tau2).^2./(1+(tau*f).^2);

mu_r = mu_p - 1j*mu_pp;
% plot(log10(f), log10(abs(mu_r)))

Ts = 0;
data = idfrd(mu_r, f, Ts);

% 1. Créer un modèle TF vide (temps continu)
np = 2;     % nombre total de pôles
nz = 2;     % nombre total de zéros (dont celui à l’origine)
Ts = 0;     % temps continu

% 2. Définir les coefficients de numérateur et dénominateur
%    (exemple : numérateur = [1 0] => zéro à l’origine)
% init_sys = idtf([1 1], [1 1]);  % structure de départ

% % 3. Geler le zéro à l'origine (on le fixe)
% init_sys.Structure.Numerator.Free = [false false];  % [z0 z1], z0 = fixe
% init_sys.Structure.Denumerator.Free = [false false];  % [z0 z1], z0 = fixe

% 4. Lancer l'identification avec structure imposée
model = tfest(data, init_sys);

% 5. Tracer la réponse
bode(data, model), grid on
legend('FRF mesurée', 'Modèle avec intégrateur')