function dydt = fv19(t, y, n, prod, degr, K, MOI, V0, convFac, ss)
%Returns differential equations at point y, no replication.

%- v13: includes basal PRM
%- v14: no lambda ODE -- lambda(t) only; AvgInj -- uses derived expression
%- v15: Dt added back in as a parameter
%- v16: Bennett prefactors
%- v17: CII degr. phenomenological, multimer degradation
%- v18: phenomenological model (total prot only, no Bennet prefac.)
%- v19: Scaling param for cII mRNA, cII tx = cro tx, PR transfer function,
%fixed minors items

if nargin < 10
    ss = 0;
end

%Initialize ODEs
dydt = zeros(6, 1);

%Correct for negative values
y(y < 0) = 0;

%Map of ODEs:
%dydt = [   
%           [cI];
%           [cro];
%           [cII];
%           [CI];
%           [Cro];
%           [CII];
%       ];

[rcI_PRM, acI_PRM, r_PRE, rCI, rcro, rCro, rcII, rCII] = deal(prod(1),...
    prod(2), prod(3), prod(4), prod(5), prod(6), prod(7), prod(8));

[kdil, kcI, kCI, kcro, kCro, kcII, kCII, kM] = deal(degr(1), degr(2), ...
    degr(3), degr(4), degr(5), degr(6), degr(7), degr(8));

[nPRM_CIu, nPRM_CId, nPRM_Cro, nPRE, nCro_Cro, nCro_CI, nCII_Cro, nCII_CI, ...
    nDeg_CII] = ...
    deal(n(1), n(2), n(3), n(4), n(5), n(6), n(7), n(8), n(9));

[KPRM_CIu, KPRM_CId, KPRM_Cro, KPRE, KCro_Cro, KCro_CI, KCII_Cro, KCII_CI, ...
    KDeg_CII] = ...
    deal(K(1), K(2), K(3), K(4), K(5), K(6), K(7), K(8), K(9));

%Viral copy number
V = V0*exp(kdil*t);
if ss ~= 1
    lambda = MOI*convFac/V;
elseif ss == 1
    lambda = MOI*convFac/V0;
end

%Prefactors----------------------------------------------------------------
PRM_prod_norm1 = 1/...
    (1 + (y(4)/KPRM_CIu)^nPRM_CIu + (y(4)/KPRM_CId)^nPRM_CId + ...
    (y(5)/KPRM_Cro)^nPRM_Cro);
PRM_prod_norm2 = (y(4)/KPRM_CIu)^nPRM_CIu/...
    (1 + (y(4)/KPRM_CIu)^nPRM_CIu + (y(4)/KPRM_CId)^nPRM_CId + ...
    (y(5)/KPRM_Cro)^nPRM_Cro);
PRE_prod_norm = (y(6)^nPRE)/(KPRE^nPRE + y(6)^nPRE);
Cro_prod_norm = 1/(1 + (y(5)/KCro_Cro)^nCro_Cro + (y(4)/KCro_CI)^nCro_CI);
CII_prod_norm = 1/(1 + (y(5)/KCII_Cro)^nCII_Cro + (y(4)/KCII_CI)^nCII_CI);
CII_deg_norm = 1/(1 + (y(6)/KDeg_CII)^nDeg_CII);
 
%[1] cI
dydt(1) = lambda*rcI_PRM*PRM_prod_norm1 + lambda*rcI_PRM*acI_PRM*PRM_prod_norm2 + ...
    lambda*r_PRE*PRE_prod_norm - ...
    kcI*y(1) - kdil*y(1);

%[2] cro
dydt(2) = lambda*rcro*Cro_prod_norm - kcro*y(2) - kdil*y(2);

%[3] cII
dydt(3) = lambda*rcII*CII_prod_norm - kcII*y(3) - kdil*y(3);

%[4] CI
dydt(4) = rCI*y(1) - kCI*y(4) - kdil*y(4);

%[5] Cro
dydt(5) = rCro*y(2) - kCro*y(5) - kdil*y(5);

%[6] CII
dydt(6) = rCII*y(3) - kCII*CII_deg_norm*y(6) - kdil*y(6);

if lambda < 0
   error('\lambda < 0!'); 
end

% if any(~isreal(dydt))
%     %keyboard
%     %dydt = real(dydt); 
%     warning('Complex ODE!');
% end;

end