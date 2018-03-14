%% function bracketing to get a bracket trial for alpha

function [retcode, ax, af, ag, alpha] = bracketing(x, f, g, d, alpha0, sigma, rho, tau1, tau2,tau3, lowerbd);

%% a                  one value of alpha bracket
%% b                  another value of alpha bracket
%% a1                 lower limit of updated alpha
%% b1                 upper limit of updated alpha
%% mu                 upper bound of initial guess for alpha
%% f_zero             alpha function at alpha = 0
%% df_zero            derivative of alpha function at alpha = 0
%% f_a                alpha function value at alpha = a
%% df_a               derivative of alpha function at alpha = a 
%% f_b                alpha function at b
%% df_b               derivative of alpha function at alpha = b 
   alpha_prev = 0.0;  %% previous value of alpha, initially 0
%% f_a_prev           alpha function value at prev alpha value
%% df_a_prev          alpha function value at prev alpha value
%% f_alpha_prev  
%% df_alpha_prev
%% f_alpha
%% df_alpha;

n = size(x,1);
f_alpha_prev = f; f_zero = f;
df_alpha_prev = d'*g; df_zero = df_alpha_prev;

mu = min(alpha0, (lowerbd - f_zero)/(rho*df_zero));
alpha = mu/2;  

k = 1;
while k < 20;

  ax = x + alpha*d;
  
  if alpha > 0.99*alpha0;
    f_alpha = 1.0e20;
  else;
    f_alpha = fun(ax);
  end;

  if f_alpha <= lowerbd; %% This really should never occur.
     retcode = -1; 
     af = f_alpha;
     ag = grad(ax);
     return;
  end;
  
  if f_alpha > f_zero + alpha*rho*df_zero | f_alpha >= f_alpha_prev;
    a = alpha_prev; f_a = f_alpha_prev; df_a = df_alpha_prev;
    b = alpha;     f_b = f_alpha;
    %%Are these gradients reasonable?
    [retcode, ax, af, ag, alpha] = sectioning(x, ax, g, d, alpha, rho, tau2, tau3, sigma, a, b, f_zero, df_zero, f_a, f_b, df_a, 0, alpha0);
    return;
  end;
   
  ag = grad(ax);
  df_alpha = d'*ag;

  if abs(df_alpha) <= -sigma*df_zero;
    %% Do we need to add this line?  
    af = f_alpha;
    retcode = 0;
    return;
  end;

  if df_alpha >= 0.0;
    a = alpha;      f_a = f_alpha;      df_a = df_alpha;
    b = alpha_prev; f_b = f_alpha_prev; df_b = df_alpha_prev;
    [retcode, ax, af, ag, alpha] = sectioning(x, ax, g, d, alpha, rho, tau2, tau3, sigma, a, b, f_zero, df_zero, f_a, f_b, df_a, 0, alpha0);
    return;
  end;
  
  if mu <= (2*alpha - alpha_prev);
    alpha_prev = alpha;
    f_alpha_prev = f_alpha;
    df_alpha_prev = df_alpha;
    alpha = mu;
  else;
    %% Store value of alpha before interpolation.
    tempval = alpha;
    %% Prepare for interpolation.
    a = alpha_prev; f_a = f_alpha_prev; df_a = df_alpha_prev;
    b = alpha;      f_b = f_alpha;      df_b = df_alpha;
    a1 = 2*alpha - alpha_prev; 
    b1 = min(mu, alpha + tau1*(alpha - alpha_prev));
    %% Do interpolation.
    alpha = interpolation(f_a, df_a, f_b, a, b, a1, b1, alpha);
     
    %% Prepare for next round.
    alpha_prev = tempval;
    f_alpha_prev = f_alpha;
    df_alpha_prev = df_alpha;
  end;
   
  k = k+1;
   
end;

retcode = -1;

