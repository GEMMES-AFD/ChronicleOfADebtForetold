#Intermediate Variables--------------------------------------------------------------------------------------------------------------------------------------------------  

##intermediate variables

#Non Financial Corporations - Intermediate --------------------------------------------------------------------------------------

gk = (ikf + FDIgreen/pk)/kf - deltaf                   #Growth rate of the net capital stock of NFCs (*Real*, *%*)
yp = ye + ivd                                          #Total production (*Real*, *%*)

ivd = betaivd*(vd - v)                                 #Desired investment in inventories (*Real*, *Flow*)
vd = alphav*ye                                         #Desired inventories (*Real*, *Stock*)

ypd = yp - im                                          #Domestic production (*Real*, *Flow*)

Yd = Con + IC + Ik + Iktr + X                          #Aggregate demand (*Nominal*, *Flow*)
yd = Con/pc + IC/pi + Ik/pk + iktr + X/px              #Aggregate demand (*Real*, *Flow*)

Con = Ch + Cg                                          #Total consumption of NFCs products (*Nominal*, *Flow*)

IC = pi*(icf + icg + icb)                              #Total intermediate consumption of NFCs products (*Nominal*, *Flow*)
icf = lambdaicf*yp                                     #NFCs' intermediate consumption demand (*Real*, *Flow*)

Ik = pk*ikf + Ikh + Ikb + Ikg + FDIgreen               #Total investment of NFCs products (*Nominal*, *Flow*)
ikfTar = (kappa0 + kappa1*(rf - pdot/p))*kf            #Target NFCs' investment demand (*Real*, *Flow*) 
kappa0 = ((1/(1 + exp(kappa01*t-kappa02)))*(kappa03 - kappa03*kappa04) +  kappa03*kappa04)


# lambdaiktr = (1/(1.0 + exp(-lambdatr0*(t - lambdatr1))))*lambdatr2
iktr = K_0*(delta_tr + alpha_tr* max(0.0000001, (t - 4))^gamma_tr * exp(- beta_tr*(t - 4)))
# lambdaiktr*K_0 + adj_iktr                    #NDC investment  (*Real*, *Flow*)
Iktr = iktr*pktr                                                                #NDC investment (*Nominal*, *Flow*)
adj_iktr_tar = (1/(1.0 + exp(-lambdatr0_adj*(t - lambdatr1_adj))))*lambdatr2_adj*max((0.013*GDP-Iktr)/pktr, 0) 

IM = im*pw*en                                                                                 #Total imports (*Nominal*, *Flow*) 
im = sigmamc*(Con/pc) + sigmamic*(IC/pi) + sigmamk*(Ik/pk) + sigmamktr*Iktr/pk                #Total imports (*Real*, *Flow*) 
sigmamcTar = sigmapcVar*(p/(pw*en*(1 + taum)))^epsilon1c + sigmaac*(aw/a)^epsilon2c           #Target propensity to import final consumption goods (*Real*, *%*) 
sigmamicTar = sigmapicVar*(p/(pw*en*(1 + taum)))^epsilon1ic + sigmaaic*(aw/a)^epsilon2ic      #Target propensity to import intermediate consumption goods (*Real*, *%*)
sigmamkTar = sigmapkVar*(p/(pw*en*(1 + taum)))^epsilon1k + sigmaak*(aw/a)^epsilon2k           #Target propensity to import investment goods (*Real*, *%*)
sigmamktrTar = sigmamktr0 + sigmaaktr*(awgr/agr)^epsilon2ktr                                     #Target propensity to import transition goods (*Real*, *%*)
taumT=max(0.064,0.064*(1+maxTaum*tanh(speedTaum*(triggerTaum-reserves))))

sigmapcVar = (1/(1 + exp(sigmamSpeed*t-sigmamInit)))*(sigmapc - sigmapc*sigmapcNew) +  sigmapc*sigmapcNew                      #Shock on the target propensity to import final consumption goods (*Real*, *%*)
sigmapicVar = (1/(1 + exp(sigmamSpeed*t-sigmamInit)))*(sigmapic - sigmapic*sigmapicNew) +  sigmapic*sigmapicNew#               #Shock on the target propensity to import intermediate consumption goods (*Real*, *%*)
sigmapkVar = (1/(1 + exp(sigmamSpeed*t-sigmamInit)))*(sigmapk - sigmapk*sigmapkNew) +  sigmapk*sigmapkNew#                     #Shock on the target propensity to import investment goods (*Real*, *%*)

X = xrO*pO*en + sigmaxn*GDPw*pwx*en                                                                     #Total exports (*Nominal*, *Flow*)
x=X/px                                                                                                  #Total exports (*Real*, *Flow*)
sigmaxnTar = sigmaxnpVar*(pwx*en/(p*(1+tauCBAM)))^epsilonxn1 + sigmaxna*((a/aw)^epsilonxn2)             #Target propensity to export non-oil and coal goods (*Real*, *%*)

sigmaxnpVar = (1/(1 + exp(sigmaxnSpeed*t-sigmaxnInit)))*(sigmaxnp - sigmaxnp*sigmaxnpNew) +  sigmaxnp*sigmaxnpNew#             #Shock on the target propensity to export non-oil and coal goods (*Real*, *%*)
tauCBAM = tauCBAM0*(ktr/kf)^tauCBAM1                                                                                          #Carbon Border Adjustment Mechanism Tax charged by the EU (*Real*, *%*)

GDP = Ch + Cg + Comh + Insh + PSg + Ik + Iktr + X - IM                         #Gross Domestic Product (*Nominal*, *Flow*)
gdp = Con/pc + Ik/pk + iktr + x - im                            #Gross Domestic Product (*Real*, *Flow*)

pd = (1 + mu)*huc                                                #Desired production price level (*Nominal*, **)
mu = mu0 - mu1*(v/ye - alphav)                                   #Mark-up over historical unitary cost (*Nominal*, *%*)
uc = ((1 + thetawf)*Wf*Lf + pi*icf + tauyf*p*ypd)/ypd            #Unitary cost (*Nominal*, **)

pc = (p*(1 - sigmamc) + pw*en*(1 + taum)*sigmamc)*(1 + tauvat + tauothc + taugreenc)                  #Price level of final consumption goods (*Nominal*, **)
pi = (p*(1 - sigmamic) + pw*en*(1 + taum)*sigmamic)*(1 + tauothi + taugreenic)                        #Price level of intermediate consumption goods (*Nominal*, **)
pk = (p*(1 - sigmamk) + pw*en*(1 + taum)*sigmamk)*(1 + tauothk)                                       #Price level of investment goods (*Nominal*, **)
px = X/(xrO + sigmaxn*GDPw)                                                                           #Price level of export goods (*Nominal*)
pktr = (p*(1 - sigmamktr) + pw*en*(1 + taumtr)*sigmamktr)*(1 + tauothktr)                             #Price level of transittion goods (*Nominal*, **)

taugreenic = shrGrTax*shrGrIC*iktr*pktr/((p*(1 - sigmamic) + pw*en*(1 + taum)*sigmamic)*(1 + tauothi)*(icf+icb+icg))  #Green tax rate levied on intermediate goods  (*%*, **)    
taugreenc = shrGrTax*(1-shrGrIC)*iktr*pktr/Con                                                              #Green tax rate levied on consumption goods   (*%*, **)

er = pw*en/p                         #Real exchange rate (*Real*, **)
pm = pw*en                           #Price level of imported goods (*Nominal*, **)

L = Lf + Lg + Lb                     #Total employment (*Real*, **)
Lf = ypd/(a)                           #Employment in NFCs (*Real*, **)
atr = atr0*(ktr/kf)^atr1            #Shock to labour intensity due to investments in transition (*Real*, **)

unem = 1 - L/LFo                     #Unemployment rate (*Real*, *%*)

GOSf = Yd - IM - Tm - tauyf*p*ypd - Tp - Tvat - pi*icf - Insf - Comf - (1 + thetawf)*Wf*Lf                   #Gross Operating Surplus of NFCs (*Nominal*, *Flow*)
GOSh = thetaGh*GOSf                  #Gross Operating Surplus redistributed from NFCs to Households (*Nominal*, *Flow*)
GOSg = thetaGg*GOSf                  #Gross Operating Surplus redistributed from NFCs to the Government (*Nominal*, *Flow*)

MIh = betaHmi*(GOSf)                 #Mixed-income redistributed from NFCs to Households (*Nominal*, *Flow*)

GFf = GOSf - Roy - ildf*Ldf - ilfxb*Lfxfb*en - ilfxfw*Lfxfw*en + idep*(1 - mdf)*Df - (MIh + GOSh + GOSg)     #Gross profits of NFCs (*Nominal*, *Flow*)
Ff = (1 - tauf)*GFf                  #Net profits of NFCs (*Nominal*, *Flow*)
FNf = Ff - Dfxfdot*en - Dfdot        #Net profits net of deposit accumulation (*Nominal*, *Flow*)

rf = Ff/(pk*kf)                      #Profit rate of NFCs (*Nominal*, *%*)

DIVf = (1 - sf)*FNf                                             #Dividends paid by NFCs (*Nominal*, *Flow*)
DIVfw = max(0,(ipsilon0w + ipsilon1w*(xrO*pO)/(GDP/en))*DIVf)          #NFCs' dividends paid to the Rest of the World (*Nominal*, *Flow*)
DIVfg = max(0,(ipsilon0g + ipsilon1g*(xrO*pO)/(GDP/en))*DIVf)          #NFCs' dividends paid to the Government (*Nominal*, *Flow*)
DIVfh = DIVf - DIVfw - DIVfg                                    #NFCs' dividends paid to Households (*Nominal*, *Flow*)

REf = sf*FNf - Othf                  #Retained earnings of NFCs (*Nominal*, *Flow*)

TFNF = pk*ikf - REf                  #Total financing needs of NFCs (*Nominal*, *Flow*)


#Financial Corporations - Intermediate -----------------------------------------------------------------------------------------------------------------------

Ybr = Ins + Com                  #Production of FCs other than FISIM (*Nominal*, *Flow*)

Com = Comh + Comf                #Total commissions paid to FCs (*Nominal*, *Flow*)
Comh = comH*Ldh                  #Commissions paid by Households (*Nominal*, *Flow*)
Comf = comF*Ldf                  #Commissions paid by NFCs (*Nominal*, *Flow*)

Ins = Insh + Insf                #Total demand for insurance services of FCs (*Nominal*, *Flow*) 
Insh = InsH*krh*pk               #Insurance services demanded by Households (*Nominal*, *Flow*)
Insf = InsF*kf*pk                #Insurance services demanded by NFCs  (*Nominal*, *Flow*)

icb = lambdaicb*Lb               #FCs' intermediate consumption demand (*Real*, *Flow*)
lambdaicb = ((1/(1 + exp(lambdaicb1*t-lambdaicb2)))*(lambdaicb3 - lambdaicb3*lambdaicb4) +  lambdaicb3*lambdaicb4)

Ikb = kappaib*Ybr                #FCs' investment demand (*Nominal*, *Flow*)

STb = (1 - fistg)*ST             #Social transfers paid by FCs (*Nominal*, *Flow*)

WSCb = (1 - phiscg)*WSC          #Workers' social contributions paid to FCs (*Nominal*, *Flow*)

GOSb = Ybr - pi*icb - tauyb*Ybr - (1 + thetawb)*Wb*Lb           #Gross operating surplus of FCs (*Nominal*, *Flow*)

GFb = ilh*Ldh + ildf*Ldf + ilfxb*Lfxfb*en + ibgdc*Bgb + Bgbtr*ibgtr - (idep*Dg + idep*(1 - mdf)*Df + idep*(1 - mdh)*Dh) - ilfxbw*Lfxbw*en + irfxb*Rfxb*en - ip*Ad - pi*icb - tauyb*Ybr - (1 + thetawb)*Wb*Lb + WSCb - STb + Ins + Com  #Gross profits of FCs (*Nominal*, *Flow*)
Fb = (1 - taub)*GFb                         #Net profits of FCs (*Nominal*, *Flow*)

OFcar = car*(Ldf + Ldh + Lfxfb*en)          #Own funds of FCs needed to accomplish the leverage regulation (*Nominal*, *Stock*)
REb = betaof*(OFcar - OFb)                  #Retained earnings of FCs (*Nominal*, *Flow*)

DIVb = Fb - REb - Othb                      #Dividends paid by FCs (*Nominal*, *Flow*)
DIVbh = DIVb                                #FCs' dividends paid to Households (*Nominal*, *Flow*)
DIVbw = 0                                   #FCs' dividends paid to the Rest of the World (*Nominal*, *Flow*)

Dd = Dg + Dh + Df                           #Total domestic currency deposits (*Nominal*, *Stock*) 

TFNB = (Ldfdot + Ldhdot + Bgbdot + Bgbtrdot) + Ikb + lr*(Dgdot + Dhdot + Dfdot) - (Dgdot + Dhdot + Dfdot + OFbdot + FDIb) - IPShdot + (Dfxbdot*en + Lfxfbdot*en - Lfxbwdot*en + Rfxbdot*en)       #Total financing needs of FCs (*Nominal*, *Flow*) 

idep = ip - md                                             #Interest rate on deposts (*Nominal*, *%*)
md = rho0 - rho1/(1 + exp(-rho2*(Ad/Dd-rho3)))             #Mark-down on Central Bank policy rate (*Nominal*, *%*)

ildf = AFC*(1 + premf)                                                               #Interest rate on domestic loans to NFCs (*Nominal*, *%*)
AFC = ((idep*Dg + idep*(1 - mdf)*Df + idep*(1 - mdh)*Dh) + ip*Ad)/(Dd + Ad)          #Average funding cost of FCs (*Nominal*, *%*) 
premfTar = zeta0 + zeta1/(1 + exp(-zeta2*((Ldf + Lfxfb*en + Lfxfw*en)/p*ypd)))       #Risk premium on loans to NFCs (*Nominal*, *%*) 

ilh = ildf*(1 + premh)                                      #Interest rate on household loans (*Nominal*, *%*) 
premhTar = chi0 + chi1/(1 + exp(-chi2*(Ldh/YDh)))           #Risk premium on household loans (*Nominal*, *%*) 

premfx = zetafx0 + zetafx1*(rsk)^zetafx2                    #Premium on FX loans with the Rest of the World (*Nominal*, *%*) 

ilfxbw = iwst + premfx                                      #Interest rate on FX loans of FCs with the Rest of the World  (*Nominal*, *%*) 
ilfxb = ilfxbw*(1 + rhofx2*premf)                           #Interest rate on FX loans of NFCs with FCs (*Nominal*, *%*) 
ilfxfw = iwst + rhofx1*premfx                               #Interest rate on FX loans with the Rest of the World (*Nominal*, *%*)

#Central Bank - Intermediate -----------------------------------------------------------------------------------------------------------------------

ipTar =   iota0 + iota1*(pdot/p - iota2)                  #Target monetary policy rate (*Nominal*, *%*)

irfx = iwst + pirfx                                       #Interest rate on FX reserves held by the Central Bank (*Nominal*, *%*)
irfxb = iwst + pirfxb                                     #Interest rate on FX reserves held by FCs (*Nominal*, *%*)

Fcb = ip*Ad + irfx*Rfxcb*en - idepcb*Dcbg                 #Central Bank profits (*Nominal*, *Flow*)

idepcb = 0.155# (*TODO AG define equation *)


#Households - Intermediate ----------------------------------------------------------------------------------------------------------------------

YDh = (1 - tauw)*WL + MIh + ESC + STg + STb - WSC - ilh*Ldh + idep*(1 - mdh)*Dh + DIVfh + DIVbh + Rem*en + GOSh - Insh - Comh + Othhh     #Households income  (*Nominal*, *Flow*)

WL = Wf*Lf + Wg*Lg + Wb*Lb                                       #Wage-bill of Households (*Nominal*, *Flow*) 
ESC = thetawf*Wf*Lf + thetawg*Wg*Lg + thetawb*Wb*Lb              #Employers' social contributions (*Nominal*, *Flow*)

WSC = ESC + phisc*WL                                             #Total workers' social contributions (*Nominal*, *Flow*)

ChTar = mpc1*YDh + mpc2*(Dh + IPSh) + Ldchdot                    #Target desired-consumption of Households (*Nominal*, *Flow*)

mpc1 = (1/(1 + exp(-lambdal0*(idep - lambdal1 - pdot/p))))*(mpcUB - mpcLB) + mpcLB       #Marginal propensity to consume out income (*Nominal*, *%*)
mpc2 = 0                                                                                 #Marginal propensity to consume out wealth (*Nominal*, *%*)   

LdchTar = thetalh*YDh                                            #Target demand for consumption loans (*Nominal*, *Flow*)

IkhTar = kappahi*YDh                                             #Target household investment demand (*Nominal*, *Flow*)
kappahi = kappah0 - kappah1*ilh - kappah2*unem                   #Households' investment to income ratio (*Nominal*, *%*)

Sh = YDh - Ch                                                    #Households savings (*Nominal*, *Flow*)

TFNH = Ikh - Sh                                                  #Total financing needs of Households (*Nominal*, *Flow*)

LFoTar = lfp*pop                           
lfp   = lfp0 + lfp1*unem
LFodot = betalf*(LFoTar - LFo)

#Government - Intermediate-----------------------------------------------------------------------------------------------------------------------

TR = Tt + Roy + GOSg + WSCg + idep*Dg + idepcb*Dcbg + DIVfg - Othg + Fcb        #Total revenue of the Government (*Nominal*, *Flow*)

Tt = Ti + Tm + Tvat + Tp + Tgr + Ty                                                   #Tax revenues of the Government (*Nominal*, *Flow*)
Ti = tauw*WL + tauf*GFf + taub*GFb                                              #Taxes on income (*Nominal*, *Flow*)
Tm = taum*IM                                                                    #Taxes on imports (*Nominal*, *Flow*)
Tvat = tauvat*Con*((1 - sigmamc)*p + sigmamc*pw*en*(1 + taum))/pc               #Value-added tax (*Nominal*, *Flow*)
Tp = tauothc*Con*((1 - sigmamc)*p + sigmamc*pw*en*(1 + taum))/pc + tauothi*IC*((1 - sigmamic)*p + sigmamic*pw*en*(1 + taum))/pi +tauothk*Ik*((1 - sigmamk)*p + sigmamk*pw*en*(1 + taum))/pk       #Other taxes on products (*Nominal*, *Flow*)                 
Tgr = taugreenc*Con*((1 - sigmamc)*p + sigmamc*pw*en*(1 + taum))/pc + taugreenic*IC*((1 - sigmamic)*p + sigmamic*pw*en*(1 + taum))/pi    #Green Tax Revenue (*Nominal*, *Flow*)
Ty = tauyf*p*ypd + tauyb*Ybr                                                    #Taxes on production (*Nominal*, *Flow*)

Roy = taur*(xrO*pO*en)                                                          #Royalties (*Nominal*, *Flow*)

WSCg = phiscg*WSC                                                               #Workers' social contributions paid to the Government (*Nominal*, *Flow*)

Gt = Gp + Gip                                                                   #Total expenditures of the Government (*Nominal*, *Flow*)
Gp = Gc + (1 + thetawg)*Wg*Lg + pi*icg + Ikg + iktr*pktr + STg                 #Primary expenditures of the Government (*Nominal*, *Flow*)
Gip = ibgdc*Bg + ibgfx*Bgfx*en + ilgfx*Lgfx*en + Lgfxtr*ilgfxtr*en  + (Bgbtr + Bgwtr)*ibgtr + (debtSwapFXBgFX*ibgfx+debtSwapFXLgFXtr*ilgfxtr+debtSwapFXLgFX*ilgfx)*(1-mdds)*(1-decds)*en               #Interest payments on public debt (*Nominal*, *Flow*)

GipNoDS = ibgdc*Bg + ibgfx*Bgfx*en + ilgfx*Lgfx*en + Lgfxtr*ilgfxtr*en  + (Bgbtr + Bgwtr)*ibgtr + (debtSwapFXBgFX*ibgfx+debtSwapFXLgFXtr*ilgfxtr+debtSwapFXLgFX*ilgfx)*en               #Interest payments on public debt (*Nominal*, *Flow*)
GipExt = ibgfx*Bgfx*en + ilgfx*Lgfx*en + Lgfxtr*ilgfxtr*en  + (debtSwapFXBgFX*ibgfx+debtSwapFXLgFXtr*ilgfxtr+debtSwapFXLgFX*ilgfx)*(1-mdds)*(1-decds)*en               #Interest payments on public debt (*Nominal*, *Flow*)

GipExtNoDS = ibgfx*Bgfx*en + ilgfx*Lgfx*en + Lgfxtr*ilgfxtr*en  + (debtSwapFXBgFX*ibgfx+debtSwapFXLgFXtr*ilgfxtr+debtSwapFXLgFX*ilgfx)*en               #Interest payments on public debt (*Nominal*, *Flow*)

Gc = PSg + Cg                                                       #Total consumption of the Government (*Nominal*, *Flow*)
PSg = (1 + thetawg)*Wg*Lg + pi*icg + deltag*pk*krg                  #Government non-market production/consumption (*Nominal*, *Flow*)
CgTar = fi2*GDP*(1-maxCg*tanh(speedCg*(triggerCg-reserves)))                                                     #Government market consumption (*Nominal*, *Flow*)

Lg = etag*pop                                                       #Total employment in the public sector (*Real*, **)

icg = lambdaicg*Lg  #lambdaicg*Lg                                                  #Government intermediate consumption demand (*Real*, *Flow*)
lambdaicg = ((1/(1 + exp(lambdaicg1*t-lambdaicg2)))*(lambdaicg3 - lambdaicg3*lambdaicg4) +  lambdaicg3*lambdaicg4)

ikgTar = kappag*krg                                                 #Target Government investment demand (*Real*, *Flow*)
IkgTar = ikgTar*pk+scenInv*GDP                                      #Target Government investment demand (*Nominal*, *Flow*)

ST = (fi3*Wf*(LFo - Lg - Lf - Lb) + fi4 * Wf*pop)                  #Notional Social transfers paid to Households (*Nominal*, *Flow*)
STg = fistg*ST*(1-maxSTg*tanh(speedSTg*(triggerSTg-reserves)))                                                      #Social transfers effectively paid by the Government (*Nominal*, *Flow*)

FD = Gt - TR - PSg                                                  #Fiscal deficit of the Government (*Nominal*, *Flow*)

TFNG = FD + Dgdot + Dcbgdot + Dfxgdot*en                            #Total financing needs of the Government (*Nominal*, *Flow*)

DgTar = fi1*Gt                                                      #Target Government deposits at FCs (*Nominal*, *Stock*)
DcbgTar = fi5*Gt                                                    #Target Government deposits at the Central Bank (*Nominal*, *Stock*)

sigmafxTar = sigmaG0 + sigmaG1*CAD                                  #Target share of public debt issuance in FX (*Nominal*, *%*) 

ibgdc = ip + premgd                                                 #Interest rate on domestic public debt (*Nominal*, *%*)
premgdTar = phi0d + phi1/(exp(-phi2*((Bg + (Bgfx+debtSwapFXBgFX*(1-decds))*en + (Lgfx+debtSwapFXLgFX*(1-decds))*en + (Lgfxtr+debtSwapFXLgFXtr*(1-decds))*en +Bgtr)/GDP)))       #Risk-premium on domestic public debt (*Nominal*, *%*)#Changed

ibgfx = iwst + rhofx3*premfx                                        #Interest rate on FX public bonds (*Nominal*, *%*)
ilgfx = (1 - rhofx4)*ibgfx                                          #Interest rate on FX public loans (*Nominal*, *%*)

ilgfxtr = (1-md_lgtr)*ilgfx                                         #Interest rate on FX Green loans
ibgtr = (1-md_bgtr)*ibgdc                                           #Interest rate on domestic green bonds


#Rest of the World - Intermediate ----------------------------------------------------------------------------------------------------------------------

CAD = -(TB + IA)/GDP                            #Current account deficit (*Nominal*, *%GDP*)
TB = X - im*pw*en                               #Trade balance of the current account (*Nominal*, *Flow*)  
IA = (Rem + irfx*Rfxcb + irfxb*Rfxb + Grants - ilgfxtr*Lgfxtr - ibgfx*Bgfx - ilgfx*Lgfx - ilfxbw*Lfxbw - ilfxfw*Lfxfw-(debtSwapFXBgFX*ibgfx+debtSwapFXLgFXtr*ilgfxtr+debtSwapFXLgFX*ilgfx)*(1-mdds)*(1-decds))*en - DIVfw - DIVbw + Othw - ibgdc*Bgw - Bgwtr*ibgtr #Income account of the current account (*Nominal*, *Flow*)

Grants = shrDon*iktr*pktr/en                    #Grants  (*Nominal*, *Flow*)

Rem = sigmaRem*GDPw*pw                          #Remittances (*Nominal*, *Flow*)

FDI = varsigmafdi*ikf*pk                        #Foreign direct investment (*Nominal*, *Flow*) 
varsigmafdi = (1/(1 + exp(varsigmafdi1*t-varsigmafdi2)))*(varsigmafdi3 - varsigmafdi3*varsigmafdi4) +  varsigmafdi3*varsigmafdi4   #Shock on total FDI (*Nominal*, *%*) 

FDIf = zetaff*FDI                               #FDI in NFCs (*Nominal*, *Flow*)
FDIgreen = shrGreenField*FDIf                   #Greenfield FDI in NFCs (*Nominal*, *Flow*)
FDInonGreen = FDIf - FDIgreen                   #Non-Greenfield FDI in NFCs (*Nominal*, *Flow*)

FDIb = (1 - zetaff)*FDI                         #FDI in FCs (*Nominal*, *Flow*)

Dfx = im*pw + ibgfx*Bgfx + ilgfx*Lgfx + ilfxbw*Lfxbw + ilfxfw*Lfxfw + ilgfxtr*Lgfxtr + DIVfw/en + DIVbw/en + Rfxbdesdot + Dfxwdot + ibgdc*Bgw/en  + Bgwtr*ibgtr/en  + (debtSwapFXBgFX*ibgfx+debtSwapFXLgFXtr*ilgfxtr+debtSwapFXLgFX*ilgfx)*(1-mdds)*(1-decds)                          #Foreign currency demand (*Nominal*, *Flow*)
Sfx = X/en + Rem + Othw/en + irfx*Rfxcb + irfxb*Rfxb + FDI/en + Bgfxdot + debtSwapFXBgFXdot + Lgfxdot + debtSwapFXLgFXdot + Lfxfwdot + Lfxbwdot + Lgfxtrdot + debtSwapFXLgFXtrdot + Grants + Bgwdot/en + Bgwtrdot/en - Rfxcbdot             #Foreign currency supply (*Nominal*, *Flow*)

NIIP = -(Rfx*en + Dfxw*en - Lfxbw*en - Lfxfw*en - (Bgfx+debtSwapFXBgFX*(1-decds))*en - (Lgfx+debtSwapFXLgFX*(1-decds))*en - (Lgfxtr+debtSwapFXLgFXtr*(1-decds))*en - Bgw -Bgwtr)/GDP                #Net International Investment Position (*Nominal*, *%GDP*)
FIP = ((Lfxfw + Lfxbw + Bgfx+debtSwapFXBgFX*(1-decds) + Lgfx+debtSwapFXLgFX*(1-decds) + Lgfxtr+debtSwapFXLgFXtr*(1-decds) - Dfxw)*en)/GDP                                              #Net FX liabilities with the Rest of the World without FDI (*Nominal*, *%GDP*) 

rsk = v1*(IM/(Rfx*en))^v2                       #Country-risk premium (*Nominal*, *%*) 

iwstTar = alphapw                               #Target short-term external interest rate (*Nominal*, *%*) 

#Indicators

perCapita=GDP/(pop*en*pw)
inflation=pDot/p
reserves=(Rfx*en)/GDP
foreignDebt=en*(Bgfx+Lgfx+Lfxfb + Lfxfw+(debtSwapFXBgFX+debtSwapFXLgFX+debtSwapFXLgFXtr)*(1-decds))/GDP
privateDebt =(Ldf + Lfxfb*en + Lfxfw*en+Ldh)/GDP
pubDebt=(Bg+(Bgfx+Lgfx+Lgfxtr+(debtSwapFXBgFX+debtSwapFXLgFX+debtSwapFXLgFXtr)*(1-decds))*en)/GDP
fiscalDef=FD/GDP
hhFrag=(ilh*Ldh)/Sh
firmsFrag=(ildf*Ldf + ilfxb*Lfxfb*en + ilfxfw*Lfxfw*en)/Ff

#Other transfers - Intermediate ---------------------------------------------------------------------------------------------------------------------------------------------------------

Othf = nuf*ypd*p                 #Other transfers of NFCs (*Nominal*, *Flow*)
Othb = nub*ypd*p                 #Other transfers of FCs (*Nominal*, *Flow)
Othg = nug*ypd*p                 #Other transfers of the Government (*Nominal*, *Flow)
Othhh = nuh*ypd*p                #Other transfers of Households (*Nominal*, *Flow)
Othw = nuw*ypd*p                 #Other transfers of the Rest of the World (*Nominal*, *Flow*)                       


#Debt Swap part
# Résolution de GiPT = Gip-iktr*pktr
ds = min(1,iktr*pktr/(en*(ibgfx*Bgfx+ilgfx*Lgfx+Lgfxtr*ilgfxtr)*(1-mdds)*(1-decds)))*dsactive

debtSwapFXLgFXTarget=ds*Lgfx
debtSwapFXLgFXtrTarget=ds*Lgfxtr
debtSwapFXBgFXTarget=ds*Bgfx

#Non-Financial Corporations - Derivatives ------------------------------------------------------------------------------------------------------------------------------------------------

yedot = betay*(yd - ye) + gk*ye                           #Change in expected sales (*Real*, *Flow*)
vdot = yp - yd                                            #Change in actual inventories (*Real*, *Flow*)

ikfdot = betaikf*(ikfTar - ikf)                           #Adjustment of NFCs' investment demand (*Real*, *Flow*)

kfdot = ikf + FDIgreen/pk - deltaf*kf                     #NFCs' capital stock accumulation (*Real*, *Flow*)
ktrdot = iktr - deltaf*ktr                              #Green capital stock accumulation (*Real*, *Flow*)

sigmamcdot = betasigmamc*(sigmamcTar - sigmamc)                   #Adjustment of the propensity to import final consumption goods (*Real*, *%*)
sigmamicdot = betasigmamic*(sigmamicTar - sigmamic)               #Adjustment of the propensity to import intermediate consumption goods (*Real*, *%*)
sigmamkdot = betasigmamk*(sigmamkTar - sigmamk)                   #Adjustment of the propensity to import investment goods (*Real*, *%*)
sigmamktrdot = betasigmamktr*(sigmamktrTar - sigmamktr)           #Adjustment of the propensity to import transition goods (*Real*, *%*)

sigmaxndot = betasigmaxn*(sigmaxnTar - sigmaxn)                   #Adjustment of the propensity to export non-oil and coal goods (*Real*, *%*)

hucdot = betahuc*(uc - huc)                               #Adjustment of the historical unitary cost (*Nominal*, **)
pdot = betap*(pd - p)                                     #Adjustment of the producer price level (*Nominal*, **)

adot = alphaa*a                                           #Domestic productivity growth (*Real*, **)
agrdot = alphaa*agr                                       #Domestic productivity growth in "green" industries (*Real*, **)

popdot = alphapop*pop                                     #Labour force growth (*Real*, **)

Wfdot = (omegaf0*(adot/a) + omegaf1*(L/LFo - omegaf2) + omegaf3*pdot/p)*Wf           #Change in wages paid by NFCs (*Nominal*, **)

Dfdot = betaDf*(etadf*Wf*Lf*(1 + thetawf) - Df)           #Adjustment of domestic NFCs' deposits (*Nominal*, *Flow*)
Dfxfdot = betaDfx*(etadfxf*(Lfxfb + Lfxfw) - Dfxf)        #Adjustment of FX NFCs' deposits (*Nominal*, *Flow*)

Lfxfbdesdot = etalfxfb*(TFNF/en)                                                 #NFCs' desired FX loans demand with FCs (*Nominal*, *Flow*)
Lfxfbdot = (1 - ratBFX)*Lfxfbdesdot                                              #NFCs' FX loans with FCs (*Nominal*, *Flow*)

Lfxfwdot = (1 - ratFFX)*etalfxfw*(TFNF/en)                                       #NFCs'FX loans with the Rest of the World (*Nominal*, *Flow*)
ratFFX = 1/(1.0 + exp(-betariskFFX*(rsk - MPFFX)))*(UBFFX - LBFFX) + LBFFX       #Credit rationing on NFCs' desired FX loans demand with the Rest of the World (*Nominal*, *%*)

Ldfdot = TFNF - Lfxfbdot*en - Lfxfwdot*en - FDInonGreen            #NFCs' domestic currency loans (*Nominal*, *Flow*) 


#Financial Corporations - Derivatives ------------------------------------------------------------------------------------------------------------------------------------------------

krbdot = Ikb/pk - deltab*krb                                 #FCs' capital stock accumulation (*Real*, *Flow*)

Lbdot = etab*Lb                                              #Change in FCs' employment demand (*Real*, **)

Wbdot = (omegab0*(adot/a) + omegab1*pdot/p)*Wb               #Change in wages paid by FCs (*Nominal*, **)

Rfxbdesdot = Lfxbwdesdot - Dfxbdot - Lfxfbdot                #FCs' desired FX reserves accumulation to accomplish with the FX no-open position condition (*Nominal,*Flow*)
Rfxbdot = Rfxdot - Rfxcbdot                                  #FCs' actual FX reserves accumulation to accomplish with the FX no-open position condition (*Nominal,*Flow*)

Dfxbdot = betadfxb*(etadbfx*Lfxbw - Dfxb)                    #Adjustment of FX FCs' deposits (*Nominal*, *Flow*)

Lfxbwdesdot = (etalxfbw*OFbdot/en + Lfxfbdesdot)                                 #FCs' desired FX loans demand with the Rest of the World (*Nominal*, *Flow*)
Lfxbwdot = (1 - ratBFX)*Lfxbwdesdot                                              #FCs'FX loans with the Rest of the World (*Nominal*, *Flow*)
ratBFX = 1/(1.0 + exp(-betariskBFX*(rsk - MPBFX)))*(UBBFX - LBBFX) + LBBFX       #Credit rationing on FCs' desired FX loans demand with the Rest of the World (*Nominal*, *%*)

Bgbdot = Bgdot - Bgwdot                                      #FCs' purchase of domestic public bonds (*Nominal*, *Flow*)
Bgbtrdot = (1-shrGrBw)*Bgtrdot                                #Domestic Green bonds purchased by FCs


Rddot = lr*((Dgdot + Dhdot + Dfdot))                         #Cash plus bank reserves retention (*Nominal*, *Flow*)

OFbdot = REb                                                 #Change in FCs' own funds (*Nominal*, *Flow*)

Addot = max(TFNB, -Ad)                                       #Liquidity advances (*Nominal*, *Flow*)

premfdot = betapremf*(premfTar - premf)                      #Adjustment of the risk premium on loans to NFCs (*Nominal*, *%*) 
premhdot = betapremh*(premhTar - premh)                      #Adjustment of the risk premium on loans to Households (*Nominal*, *%*)

#Central Bank - Derivatives ------------------------------------------------------------------------------------------------------------------------------------------------

# Rfxcbdot = max(0,(sigmaRfxb0+sigmaRfxb1*IM/(Rfx*en))*(sigmaRfxb*im*pw - Rfxcb))                 #Central Bank's intervention in the FX market (*Nominal*, *Flow*)
Rfxcbdot = max(sigmaRfxb*im*pw + Rfxcb*shareFX*(1+tanh(speedFX*(triggerFX-reserves))) - Rfxcb, 0)                     #Central Bank's intervention in the FX market (*Nominal*, *Flow*)
ipdot = betaip*(ipTar - ip)                                    #Adjustment of the monetary policy rate (*Nominal*, *%*)


#Households - Derivatives ------------------------------------------------------------------------------------------------------------------------------------------------

Chdot = betacon*(ChTar - Ch)                               #Convergence of Households' consumption (*Nominal*, *Flow*)

Ikhdot = betaIh*(IkhTar - Ikh)                             #Adjustment of the Households' investment demand (*Real*, *Flow*)
krhdot = Ikh/pk - deltah*krh                               #Households' capital stock accumulation (*Real*, *Flow*)

Ldhdot = Ldchdot + Ldihdot                                 #Total households' loans demand (*Nominal*, *Flow*)
Ldihdot = thetal3*Ikh                                      #Mortgage loans (*Nominal*, *Flow*)
Ldchdot = betaLdch*(LdchTar - Ldch)                        #Adjustment of consumption loans demand (*Nominal*, *Flow*)
thetalhdot = 0                                             #Adjustment of consumption loans demand to households income ratio (*Nominal*, *Flow*)

IPShdot = zetaitr*WL                                       #Insurance, pensions and standarised guaranteed schemes accumulation (*Nominal*, *Flow*)

Dhdot = Sh - Ikh + Ldhdot - IPShdot                        #Households' accumulation of domestic deposits (*Nominal*, *Flow*)

#Government - Derivatives ------------------------------------------------------------------------------------------------------------------------------------------------

Cgdot = betaCg*(CgTar - Cg)                                  #Adjustment of the Government market-consumption (*Nominal*, *Flow*)

Ikgdot = betaIkg*(IkgTar - Ikg)                              #Adjustment of the Government investment demand (*Real*, *Flow*)
krgdot = Ikg/pk - deltag*krg                                 #Governments' capital stock accumulation (*Real*, *Flow*)

Wgdot = ((omegag0*(adot/a) + omegag1*pdot/p)*(1-maxWg*tanh(speedWg*(triggerFX-reserves))))*Wg               #Change in wages paid by the Government (*Nominal*, *Flow*)

Dgdot = betaDg*(DgTar - Dg)                                  #Government accumulation of domestic deposits at FCs (*Nominal*, *Flow*)
Dcbgdot = betaDcbg*(DcbgTar - Dcbg)                          #Government accumulation of domestic deposits at the Central Bank (*Nominal*, *Flow*)
Dfxgdot = betaDfxg*(etadfxg*(Lgfx +debtSwapFXLgFX+ Bgfx +debtSwapFXBgFX) - Dfxg)            #Government accumulation of FX deposits (*Nominal*, *Flow*)          #Government accumulation of FX deposits (*Nominal*, *Flow*)

sigmafxdot = betasigmafx*(sigmafxTar-sigmafx)                #Adjustment of the share of public debt issuance in FX (*Nominal*, *%*) 

Bgfxdot = -zetabgfx*(sigmafx*TB/en) - debtSwapFXBgFXdot                              #Government FX bonds issuance (*Nominal*, *Flow*)

Lgfxdot = -(1 - zetabgfx)*(sigmafx*TB/en)-debtSwapFXLgFXdot  #Government FX loans (*Nominal*, *Flow*)
# Lgfxdot = -(1 - zetabgfx)*(sigmafx*TB/en)                  #Government FX loans (*Nominal*, *Flow*)

Lgfxtrdot = shrGrL*shrGrLFx*iktr*pktr/en-debtSwapFXLgFXtrdot                         #Government FX green loans (*Nominal*, *Flow*)

Bgtrdot = shrGrL*(1-shrGrLFx)*iktr*pktr                       #Government domestic green bonds (*Nominal*, *Flow*)

Bgdot = TFNG - Bgfxdot*en - Lgfxdot*en - Lgfxtrdot*en - debtSwapFXLgFXdot*en - debtSwapFXLgFXtrdot*en - debtSwapFXBgFXdot*en -Bgtrdot - shrDon*iktr*pktr/en       #Government domestic bonds issuances (*Nominal*, *Flow*)      #Government domestic bonds issuances (*Nominal*, *Flow*)

premgddot = betapremgd*(premgdTar - premgd)                  #Adjustes of the risk-premium on domestic public debt (*Nominal*, *%*)

#Debt swap part
debtSwapFXLgFXdot=betadebtSwapFXLgFX*(debtSwapFXLgFXTarget-debtSwapFXLgFX)
debtSwapFXLgFXtrdot=betadebtSwapFXLgFXtr*(debtSwapFXLgFXtrTarget-debtSwapFXLgFXtr)
debtSwapFXBgFXdot=betadebtSwapFXBgFX*(debtSwapFXBgFXTarget-debtSwapFXBgFX)



#Rest of the World - Derivatives ------------------------------------------------------------------------------------------------------------------------------------------------

Bgwdot = -zetabg*TB                                     #Rest of the World's purchase of domestic public bonds (*Nominal*, *Flow*)

Bgwtrdot = shrGrBw*Bgtrdot                                    #Green bonds purchased by the RoW

Rfxdot = (TB/en + IA/en + FDI/en + Bgfxdot + Bgwdot/en + Bgwtrdot/en+ Lgfxdot + Lfxfwdot + Lfxbwdot + Lgfxtrdot +debtSwapFXLgFXdot+debtSwapFXBgFXdot+debtSwapFXLgFXtrdot + Grants - Dfxwdot)    #Total FX reserves accumulation (*Nominal*, *Flow*) 

Dfxwdot = Dfxfdot + Dfxbdot + Dfxgdot                   #Total FX deposits accumulation (*Nominal*, *Flow*)

endot = betaen*((Dfx - Sfx)/Sfx)                        #Change in nominal exchange rate (*Nominal*, **)

iwstdot = betaiwst*(iwstTar - iwst)                     #Adjustment in short-term external interest rate (*Nominal*, *%*)

awdot = alphaw*aw                                       #Foreign productivity growth rate (*Real*, **)
awgrdot = alphaw*awgr                                   #Foreign productivity growth rate in "green" industries (*Real*, **)

pwdot = alphapw*pw                                      #Foreign inflation rate (*Nominal*, **)
pwxdot = alphapw*pwx                                    #Foreign GDP deflator inflation (*Nominal*, **)   
pOdot = alphapO*pO                                      #Foreign inflation in oil and coal prices (*Nominal*, **)

GDPwdot = alphagw*GDPw                                  #Foregin GDP growth rate (*Real*, **)

#Debt Swap part
DSFXdot = debtSwapFXLgFX

#Time Derivatives------------------------------------------------------------------------------------------------------------------------------------------------------------------

##time derivatives

ye = yedot                               #Real expected growth (*Real*, **)
v = vdot                                 #Inventories accumulation (*Real*, *Flow*)
kf = kfdot                               #NFCs' capital stock accumulation (*Real*, *Flow*)
ktr = ktrdot                           #"Green" capital stock accumulation (*Real*, *Flow*)
ikf = ikfdot                             #Change in NFCs' investment demand (*Real*, **)
sigmamc = sigmamcdot                     #Change in propensity to import final consumption goods (*Real*, *%*)
sigmamic = sigmamicdot                   #Change in propensity to import intermediate consumption goods (*Real*, *%*)
sigmamk = sigmamkdot                     #Change in propensity to import investment goods (*Real*, *%*)
sigmamktr = sigmamktrdot                 #Change in propensity to import transition goods (*Real*, *%*)
sigmaxn = sigmaxndot                     #Change in propensity to export non-oil and coal (*Real*, *%*) 
huc = hucdot                             #Change in the historial unitary cost (*Nominal*, **)
p = pdot                                 #Producer price inflation (*Nominal*, **)
Wf = Wfdot                               #Growth of the average wage paid by NFCs (*Nominal*, **)
Df = Dfdot                               #NFCs' domestic deposits accumulation (*Nominal*, *Flow*)
Dfxf = Dfxfdot                           #NFCs' FX deposits accumulation (*Nominal*, *Flow*)
Lfxfb = Lfxfbdot                         #Change in NFCs' FX loans with FCs (*Nominal*, *Flow*)
Lfxfw = Lfxfwdot                         #Change in NFCs' FX loans with the Rest of the World (*Nominal*, *Flow*)
Ldf = Ldfdot                             #Change in NFCs' domestic currency loans (*Nominal*, *Flow*)
Lb = Lbdot                               #Change in FCs' employment level (*Real*, **)
krb = krbdot                             #FCs' capital stock accumulation (*Real*, *Flow*)
Wb = Wbdot                               #Growth in the average wage paid by FCs (*Nominal*, **)
Rfxb = Rfxbdot                           #FCs' FX reserves accumulation (*Nominal*, *Flow*)
Dfxb = Dfxbdot                           #FCs' FX deposits accumulation (*Nominal*, *Flow*)
Lfxbw = Lfxbwdot                         #Change in FCs' FX loans with the Rest of the World (*Nominal*, *Flow*)
Ad = Addot                               #Change in liquidity advances (*Nominal*, *Flow*)
Bgb = Bgbdot                             #Government domestic currency bonds purchased by FCs (*Nominal*, *Flow*)
Rd = Rddot                               #Cash plus bank reserves accumulated at the Central Bank (*Nominal*, *Flow*)
OFb = OFbdot                             #FCs' own funds accumulation (*Nominal*, *Flow*)
premf = premfdot                         #Change in NFCs' premium on domestic currency loans rate (*Nominal*, *%*)
premh = premhdot                         #Change in Households' premium on domestic currency loans rate (*Nominal*, *%*)
ip = ipdot                               #Change in monetary policy rate (*Nominal*, *%*)
Ch = Chdot                               #Change in desired Households' consumption (*Nominal*, **)
krh = krhdot                             #Households' capital stock accumulation (*Real*, *Flow*)
Ikh = Ikhdot                             #Change in Households' investment demand (*Nominal*, **)
Dh = Dhdot                               #Households' domestic deposits accumulation (*Nominal*, *Flow*)
Ldh = Ldhdot                             #Change in Households' domestic currency loans (*Nominal*, *Flow*)
IPSh = IPShdot                           #Households' insurance, pensions and SGS accumulation (*Nominal*, *Flow*)
Ldih = Ldihdot                           #Change in Households' mortage loans (*Nominal*, *Flow*)
thetalh = thetalhdot                     #Change in consumption loans to Households' disposable income ratio (*Nominal*, *%*)
Wg = Wgdot                               #Growth in the average wage paid by the Government (*Nominal*, **)
Ikg = Ikgdot                             #Change in Government's investment demand (*Nominal*, **)
Cg = Cgdot                               #Change in Government's market-consumption (*Nominal*, **)
krg = krgdot                             #Government's capital stock accumulation (*Real*, *Flow*)
Dg = Dgdot                               #Government's domestic currency deposits accumulation at FCs (*Nominal*, *Flow*)
Dcbg = Dcbgdot                           #Government's domestic currency deposits accumulation at the Central Bank (*Nominal*, *Flow*)
Dfxg = Dfxgdot                           #Government's FX deposits accumulation (*Nominal*, *Flow*)
Bg = Bgdot                               #Government domestic currency bonds issuance (*Nominal*, *Flow*)
Bgfx = Bgfxdot                           #Government FX bonds issuance (*Nominal*, *Flow*)
Bgw = Bgwdot                             #Government domestic currency bonds purchased by the Rest of the World (*Nominal*, *Flow*)
Lgfx = Lgfxdot                           #Change in Government's FX loans (*Nominal*, *Flow*)
Lgfxtr = Lgfxtrdot                       #Change in Green FX loans (*Nominal*, *Flow*)
Bgtr = Bgtrdot                           #Domestic Green bonds issued by the Government (*Nominal*, *Flow*)
Bgwtr = Bgwtrdot                         #Domestic Green bonds purchased by the RoW (*Nominal*, *Flow*)       
Bgbtr = Bgbtrdot                         #Domestic Green bonds purchased by FCs (*Nominal*, *Flow*)   
premgd = premgddot                       #Change in the premium on public domestic bonds rate (*Nominal*, *%*)
Rfx = Rfxdot                             #Total FX reserves accumulation (*Nominal*, *Flow*)
Dfxw = Dfxwdot                           #Total FX deposits accumulation (*Nominal*, *Flow*) 
Rfxcb = Rfxcbdot                         #Central Bank's FX reserves accumulation (*Nominal*, *Glow*)
en = endot                               #Change in the nominal exchange rate (*Nominal*, **)
a = adot                                 #Domestic labour productivity growth (*Real*, **)
aw = awdot                               #Foreign labour productivity growth (*Real*, **)
agr = agrdot                             #Domestic labour productivity growth in "green" industries (*Real*, **)
awgr = awgrdot                           #Foreign labour productivity growth in "green" industries (*Real*, **)
pw = pwdot                               #Foreign imports price level inflation (*Nominal*, **)
pwx = pwxdot                             #Foreign GDP deflator inflation (*Nominal*, **)
pO = pOdot                               #Implied oil and coal price level inflation (*Nominal*, **)
GDPw = GDPwdot                           #Foreign real GDP growth (*Real*, **)
pop = popdot                             #Labour force growth (*Real*, **)
LFo = LFodot
Ldch = Ldchdot                           #Change in Households' consumption loans (*Nominal*, *Flow*)
iwst = iwstdot                           #Change in the short-term external interest rate (*Nominal*, *%*)
sigmafx = sigmafxdot                     #Change in the share of public debt issuance in FX (*Nominal*, *%*)
xrO=-reducXrO*xrO                        #Fossil fuel exports growth (*Real*, **)
adj_iktr = 2*(adj_iktr_tar-adj_iktr)
taum=betataum*(taumT-taum)

#Debt Swap Parameters
debtSwapFXLgFX=betadebtSwapFXLgFX*(debtSwapFXLgFXTarget-debtSwapFXLgFX)
debtSwapFXLgFXtr=betadebtSwapFXLgFXtr*(debtSwapFXLgFXtrTarget-debtSwapFXLgFXtr)
debtSwapFXBgFX=betadebtSwapFXBgFX*(debtSwapFXBgFXTarget-debtSwapFXBgFX)

#Initial Values----------------------------------------------------------------------------------------------------------------------------

##initial values

ye=1674.564                           #Real expected sales (*Real*, *Flow*)
v=127.1373                            #Inventories  (*Real*, *Stock*)
kf=2213.6793921781                    #NFCs' capital stock (*Real*, *Stock*)
ktr=0                                 #Green capital stock (*Real*, *Stock*)
ikf=104.167738781951                  #NFCs' investment demand (*Real*, *Flow*)
sigmamc=0.11989                       #Propensity to import final consumption goods (*Real*, *%*)
sigmamic=0.09466                      #Propensity to import intermediate consumption goods (*Real*, *%*)                     
sigmamk=0.29089                       #Propensity to import investment goods (*Real*, *%*) 
sigmamktr=0.2870677*2                         #Propensity to import transition goods (*Real*, *%*) 
sigmaxn=0.001288757                   #Propensity to export non-oil and coal (*Real*, *%*) 
huc=0.714649985328843                 #Historial unitary cost (*Nominal*, **)
p=1.087                               #Producer price level (*Nominal*, **)
Wf=11.443                             #Average wage paid by NFCs (*Nominal*, **)
Df=82.931                             #NFCs' domestic deposits (*Nominal*, *Stock*)
Dfxf=16.73778                         #NFCs' FX deposits (*Nominal*, *Stock*)
Lfxfb=21.5441                         #NFCs' FX loans with FCs (*Nominal*, *Stock*)
Lfxfw=97.22792                        #NFCs' FX loans with the Rest of the World (*Nominal*, *Stock*)
Ldf=234.187                           #NFCs' domestic currency loans (*Nominal*, *Stock*)
Lb=0.348                              #FCs' employment level (*Real*, **)
krb=9.917121116                       #FCs' capital stock (*Real*, *Stock*)  
Wb=47.09                              #Average wage paid by FCs (*Nominal*, **)
Rfxb=40.27344                         #FCs' FX reserves (*Nominal*, *Stock*)
Dfxb=11.71839                         #FCs' FX deposits (*Nominal*, *Stock*)
Lfxbw=69.91541                        #FCs' FX loans with the Rest of the World (*Nominal*, *Stock*)
Ad=8.55                               #Liquidity advances (*Nominal*, *Stock*)
Bgb=254.865                           #Government domestic currency bonds held by FCs (*Nominal*, *Stock*)
Rd=102.482                            #Cash plus bank reserves at the Central Bank (*Nominal*, *Stock*)
OFb=105.95                            #FCs' own funds (*Nominal*, *Stock*)
premf=1.65962                         #NFCs' premium on domestic currency loans rate (*Nominal*, *%*)
premh=0.44049                         #Households' premium on domestic currency loans rate (*Nominal*, *%*)
ip=0.058                              #Monetary policy rate (*Nominal*, *%*)
Ch=686.6039807                       #Desired Households' consumption (*Nominal*, *Flow*)
krh=368.88722                         #Households' capital stock (*Real*, *Stock*)
Ikh=49.722                            #Households' investment demand (*Nominal*, *Flow*)
Dh=173.869                            #Households' domestic deposits (*Nominal*, *Stock*)
Ldh=198.497                           #Households' domestic currency loans (*Nominal*, *Stock*)
IPSh=506.068                          #Households' insurance, pensions and SGS (*Nominal*, *Stock*)
Ldih=62.111                           #Households' mortage loans (*Nominal*, *Stock*)
thetalh=0.321717066621123             #Consumption loans to Households' disposable income ratio (*Nominal*, *%*)
Wg=29.596                             #Average wage paid by the Government (*Nominal*, **)
Ikg=36.683                            #Government's investment demand (*Nominal*, *Flow*)
Cg=21.601                             #Government's market-consumption (*Nominal*, *Flow*)
krg=498.486630008789                  #Government's capital stock (*Real*, *Stock*)
Dg=50.705                             #Government's domestic currency deposits at FCs (*Nominal*, *Stock*)
Dcbg=7.943                            #Government's domestic currency deposits at the Central Bank (*Nominal*, *Stock*)
Dfxg=8.599901                         #Government's FX deposits (*Nominal*, *Stock*)
Bg=338.851                            #Government domestic currency bonds (*Nominal*, *Stock*)
Bgfx=147.8425                         #Government FX bonds (*Nominal*, *Stock*)
Bgw=83.986                            #Government domestic currency bonds held by the Rest of the World (*Nominal*, *Stock*)
Lgfx=59.73074                         #Government FX loans (*Nominal*, *Stock*)
Lgfxtr=0                              #FX Green loans (*Nominal*, *Stock*)
Bgtr = 0                              #Domestic Green bonds (*Nominal*, *Stock*)
Bgwtr = 0                             #Domestic Green bonds held by the Row  (*Nominal*, *Stock*)
Bgbtr = 0                             #Domestic Green bonds held by FCs (*Nominal*, *Stock*)
premgd=0.012                          #Premium on public domestic bonds rate (*Nominal*, *%*)
Rfx=171.8851                          #Total FX reserves (*Nominal*, *Stock*)
Dfxw=32.20152                         #Total FX deposits (*Nominal*, *Stock*)
Rfxcb=131.6116                        #Central Bank's FX reserves (*Nominal*, *Stock*)
en=1.19513                            #Nominal exchange rate (*Nominal*, **)
a=74.45                               #Domestic labour productivity (*Real*, **)
aw=80                                 #Foreign labour productivity (*Real*, **)
agr=74.45                             #Domestic labour productivity in "green" industries (*Real*, **)
awgr=80                               #Foreign labour productivity in "green" industries (*Real*, **)
pw=0.9513                             #Foreign imports price level (*Nominal*, **)
pwx=0.555204468                       #Foreign GDP deflator (*Nominal*, **)
pO=4.038202828                        #Implied oil and coal price level (*Nominal*, **)
GDPw=88952.63574                      #Foreign real GDP (*Real*, *Flow*)
pop= 49.3191                          #Scaled due to change in population definition 
LFo = 24.404942                        #Labour force (*Real*, **)
Ldch=236.386                          #Households' consumption loans (*Nominal*, *Stock*)
iwst=0.023175727                      #Short-term external interest rate (*Nominal*, *%*)
sigmafx=0.04                          #Share of public debt issuance in FX (*Nominal*, *%*)
xrO = 19.09267727                     #Fossil fuel exports (*Real*, *Flow*)
adj_iktr = 0
taum=0.064085801 
#Parameters---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

debtSwapFXLgFX=0
debtSwapFXLgFXtr=0
debtSwapFXBgFX=0

##parameters

#Non-Financial Corporations - Parameters-----------------------------------------------------------------------------------------------------------------------------------------------

betay=3                           #Speed of convergence of expected demand

alphav=0.07831128                 #Ratio of desired inventories to expected sales
betaivd=0.16300412                #Speed of convergence for desired inventory accumulation

lambdaicf=0.4063                  #Intermediate consumption technical coefficient from NFCs

#kappa0=0.0446637659861748         #Autonomous parameter in NFCs investment function

kappa01 = 0.2
kappa02 = 1
kappa03 = 0.0370265
kappa04 = 1

kappa1=0.5                        #Sensitivity of NFCs investment to real profit rate
betaikf=1                         #Speed of convergence for NFCs investment
deltaf=0.04                       #Depreciation rate of NFCs' capital stock

sigmapc=0.129665973681575         #Linear term in price effect on propensity to import final consumption goods
sigmapcNew=1                      #Size of the shock the propensity to import final consumption goods
epsilon1c=0.75425                 #Price elasticity for propensity to import final consumption goods
sigmaac=0.00025                   #Linear term in productivity effect on propensity to import final consumption goods
epsilon2c=1.5758                  #Productivity elasticity for propensity to import final consumption goods
betasigmamc=0.58575193478         #Speed of adjustment for propensity to import final consumption goods

sigmapic=0.100863727514328        #Linear term in price effect on propensity to import intermediate consumption goods
sigmapicNew=1                     #Size of the shock the propensity to import intermediate consumption goods
epsilon1ic=0.6912                 #Price elasticity for propensity to import intermediate consumption goods
sigmaaic=0.00085                  #Linear term in productivity effect on propensity to import intermediate consumption goods
epsilon2ic=2.0637                 #Productivity elasticity for propensity to import intermediate consumption goods
betasigmamic=0.64114086974        #Speed of adjustment for propensity to import intermediate consumption goods

sigmapk=0.304856100732353         #Linear term in price effect on propensity to import investment goods
sigmapkNew=1                      #Size of the shock the propensity to import investment goods
epsilon1k=0.4454                  #Price elasticity for propensity to import investment goods
sigmaak=0.00022                   #Linear term in productivity effect on propensity to import investment goods
epsilon2k=0.241                   #Productivity elasticity for propensity to import investment goods
betasigmamk=1.77651522362         #Speed of adjustment for propensity to import investment goods

sigmamSpeed=0.3                   #Speed of the shock on propensities to import
sigmamInit=4                      #Initial period of the shock on propensities to import

sigmaxnp=0.001428421              #Linear term in price effect on the propensity to export
epsilonxn1=0.6                    #Price elasticity for propensity to export
sigmaxna=0.00025                  #Linear term in productivity effect on the propensity to export
epsilonxn2=1.3745340223           #Productivity elasticity for propensity to export
betasigmaxn=1                     #Speed of adjustment for propensity to export

sigmaxnpNew=1                     #Size of the shock on the propensity to export
sigmaxnSpeed=0.4                    #Speed of the shock on the propensity to export
sigmaxnInit=5                     #Initial period of the shock on propensities to export

mu0=0.59143893                    #Autonomous parameter in mark-up function
mu1=0.01321461                    #Sensitivity of the mark-up to inventory accumulation
betahuc=15                        #Speed of convergence of historical unitary cost
betap=0.75                        #Speed of convergence of prices

thetaGh=0.09721126                #Share of NFCs' gross operating surplus distributed to Households
thetaGg=0.005745633               #Share of NFCs' gross operating surplus distributed to the Government
betaHmi=0.380038148084672         #Share of NFCs' gross operating surplus distributed as mixed income to Households

omegaf0=1                         #Sensitivity of NFCs’ wage curve to productivity
omegaf1=0                         #Sensitivity of NFCs’ wage curve to unemployment
omegaf2=0.8797                    #Employment rate as reference for wage curve
omegaf3=1                         #Sensitivity of NFCs’ wage curve to prices
thetawf=0.1642046                 #Share of NFCs' wage bill paid as social security contributions to Households

alphaa=0.02                       #Productivity growth rate
alphapop=0.005                     #Labour force growth rate

ipsilon0w=0.02             #Autonomous term on NFCs' dividends distribution towards the Rest of the World
ipsilon1w=1.095614749                       #Sensitivity of dividends distribution to fossil fuel exports towards the Rest of the World
ipsilon0g=0.019036217             #Autonomous term on NFCs' dividends distribution towards the Government
ipsilon1g=0.848266913                       #Sensitivity of dividends distribution to fossil fuel exports towards the Government

sf=0.460959487668661              #Saving rate of NFCs

etadf=0.35                        #Share of NFCs' wage bill used to determine NFCs’ target domestic deposits
betaDf=1                          #Speed of convergence of NFCs’ domestic deposits

etadfxf=0.14                      #Share of NFCs' FX loans used to determine NFCs' target FX deposits
betaDfx=1                         #Speed of convergence of NFCs' FX deposits

etalfxfb=0.04                     #Share of NFCs’ financing needs asked in the form of FX loans from FCs
etalfxfw=0.13                     #Share of NFCs’ financing needs asked in the form of FX loans from the Rest of the World

betariskFFX=30                    #Elasticity of NFCs’ rationing to country risk
UBFFX=1.11879                     #Upper bound for rationing of NFCs' FX loans
LBFFX=0.029903                    #Lower bound for rationing of NFCs' FX loans
MPFFX=0.08                        #Country risk reference for NFCs rationing

#ratFFX=0.15                       #Credit rationing on NFCs' desired FX loans demand with the Rest of the World


#Financial Corporations - Parameters-----------------------------------------------------------------------------------------------------------------------------------------------

comH=0.09065261                    #Commission generation from Households debt
comF=0.1365994                     #Commission generation from NFCs debt

InsH=0                             #Insurance generation from  Households capital stock
InsF=0                             #Insurance generation from NFCs capital stock

#lambdaicb=75.1592                  #Intermediate consumption technical coefficient from FCs
lambdaicb1 = 0.2
lambdaicb2 = 4
lambdaicb3 = 75.1592
lambdaicb4 = 1.6


etab=0.01                          #Growth rate of FCs' employment demand

kappaib=0.06862745                 #Propensity to invest out of FCs production
deltab=0.0448                      #Depreciation rate of FCs capital stock

omegab0=1                          #Sensitivity of FCs’ wage curve to productivity
omegab1=1                          #Sensitivity of FCs’ wage curve to prices
thetawb=0.1884146                  #Share of FCs' wage bill paid as social security contributions to Households

car=0.2742203                      #Capital adequacy ratio
betaof=1                           #Speed of convergence for capital adequacy ratio

etadbfx=0.0535                     #Share of FCs’ FX loans to determine FCs' target FX deposits
betadfxb=0.8                       #Speed of convergence for FCs' FX deposits

rho0=0.530451                      #Autonomous term for mark-down on policy rate for deposit rate
rho1=1                             #Linear term for mark-down on policy rate for deposit rate
rho2=1                             #Elasticity of mark-down on policy rate for deposit rate on advances over deposits                   
rho3=0.01                          #Target for Advances over Deposits for mark-down on policy rate for deposit rate

mdf=0.282092928                    #Mark-down on deposit interest rate for NFCs
mdh=0.314469636                    #Mark-down on deposit interest rate for Households

zeta0=0.915406                     #Autonomous term for risk premium of NFCs
zeta1=1                            #Linear term for risk premium of NFCs
zeta2=4                            #Elasticity of risk premium of NFCs on total debt over domestic production
betapremf=0.8507                   #Speed of convergence of NFCs premium 

chi0=0.155                         #Autonomous term for risk premium of Households
chi1=0.5                             #Linear term for risk premium of Households
chi2=1                             #Elasticity of risk premium of Households on total debt over disposable income
betapremh=1.2339                   #Speed of convergence of Households premium 

zetafx2=1                          #Autonomous term for FX risk premium
zetafx0=0.003                      #Linear term for FX risk premium
zetafx1=0.700904682321678          #Elasticity of FX risk premium on country risk

rhofx1=6.4                         #Mark-up on NFCs’ FX loans rate with the Rest of the World
rhofx2=0.185277                    #Multiplier of NFCs' FX loans rate with FCs on their domestic risk premium

# betadfxb=0.8                       #Speed of convergence for FCs' FX deposits

etalxfbw=0.2                       #Share of FCs’ own funds used for FX loans determination

betariskBFX=30                     #Elasticity of FCs’ rationing to country risk
UBBFX=1.07687                      #Upper bound for rationing of FCs in FX
LBBFX=0.37229                      #Lower bound for rationing of FCs in FX
MPBFX=0.08                         #Country risk reference for FCs rationing in FX

lr=0.1585843                       #Liquidity retention ratio

#Central Bank - Parameters-----------------------------------------------------------------------------------------------------------------------------------------------

iota0=0.04852755                   #Autonomous term for monetary policy rate     
iota1=2                            #Linear term for monetary policy rate on inflation
iota2=0.03                         #Target inflation for policy rate
betaip=0.5                         #Speed of convergence for monetary policy rate

pirfx=-0.003146882                 #Mark-up on official reserves return rate              
pirfxb=0.01935873                  #Mark-up on FCs reserve return rate

sigmaRfxb=0.74                     #Share of imports to determine FX reserves target of the Central Bank

#Households - Parameters-----------------------------------------------------------------------------------------------------------------------------------------------

phisc=0.01585536                   #Social contribution rate on households wage bill

betacon=1                          #Speed of convergence for Households consumption
lambdal0=-6                        #Elasticity of marginal propensity to consume to the real interest rate
lambdal1=0.013                     #Target real deposit rate in propensity to consume
mpcUB=0.97                         #Upper bound of the marginal propensity to consume
mpcLB=0.87                         #Lower bound of the marginal propensity to consume

betaLdch=12                        #Speed of convergence for consumption credit demand
thetal3=0.126704476891517          #Share of Households’ investment borrowed

kappah0=0.0727037178355884         #Autonomous term in Households' propensity to invest
kappah1=0                          #Sensitivity to the interest rate of Households' propensity to invest
kappah2=0                          #Sensitivity to the unemployment rate of Households' propensity to invest

betaIh=1                           #Speed of convergence for Households’ investment
deltah=3.2e-06                     #Depreciation rate of Households' capital stock

zetaitr=0.03468189                 #Share of Households' wage bill saved as insurance, pensions and SGS

lfp0 = 0.5643876
lfp1 = -0.6335344
betalf = 1

#Government - Parameters-----------------------------------------------------------------------------------------------------------------------------------------------

tauf=0.1897134                   #Tax rate on NFCs' profits                    
taub=0.1246998                     #Tax rate on FCs' profits
tauw=0.09166008                   #Tax rate on Households' wage bill

betataum= 15
tauvat= 0.1111355             #Value-added tax rate
tauothc=0                          #Other consumption tax rate 
tauothi=0.024434906                #Other intermediate consumption tax rate
tauothk=0.024434906                #Other investment tax rate

tauyf=0.0172865778595673           #Tax rate on NFCs' production
tauyb=0.03                         #Tax rate on FCs' production

taur=0.102013522823697             #Propensity of fossil fuel export to generate royalties

phiscg=1                           #Share of workers' social contributions paid to the Government

omegag0=1                          #Sensitivity of Goverment's wage curve to productivity
omegag1=1                          #Sensitivity of Goverment's wage curve to prices
thetawg=0.2721489                  #Share of Government's wage bill paid as social security contributions to Households

etag=0.04306648            #Share of population employed as public employees

#lambdaicg=23.17883                 #Intermediate consumption technical coefficient from the Government
lambdaicg1 = 0.2
lambdaicg2 = 4
lambdaicg3 = 23.17884
lambdaicg4 = 1.5


kappag=0.07                        #Government investment rate
deltag=0.035                       #Depreciation rate of the Government's capital stock
betaIkg=1                          #Speed of convergence of Government investment

fi2=0.0225                         #Share of GDP for Government market consumption
betaCg=1                           #Speed of convergence for Government market consumption

fi3= 0.45        #Share of NFCs' wages paid as unemployment benefit
fi4= 0.1682448            #0.511678    #0.5936                         #Share of NFCs' wages paid as social transfer
fistg=0.929015               #0.929015        #Share of social transfers paid by the Government

fi1=0.11209336702978144            #Share of Government expenditure to determine domestic currency deposits at FCs
betaDg=0.965                       #Speed of convergence of Government deposits at FCs

fi5=0.024740401610319002           #Share of government expenditure to determine domestic currency deposits at the central bank
betaDcbg=1                         #Speed of convergence of Government deposits at the Central Bank

etadfxg=0.0265                     #Share of FX public debt to determine Government FX deposits
betaDfxg=1                         #Speed of convergence of Government FX deposits

sigmaG0=0.145                       #Autonomous term of share of the trade balance to determine FX public debt
sigmaG1=1.5                          #Sensitivity of share of trade balance to determine FX public from the current Account deficit as a share of GDP
betasigmafx=3                      #Speed of convergence for public debt issuance in FX

zetabgfx=0.6                       #Share of FX bonds in total FX debt issuance

phi0d=0.006                        #Autonomous term of the premium of public bond rate over policy rate
phi1=0.00055                             #Linear term of the premium of public bond rate over policy rat
phi2=4                             #Elasticity of the premium of public bond rate from public debt to GDP ratio
betapremgd=1                       #Speed of convergence for domestic public debt premium

rhofx3=0.887389870407265           #Multiplier of country risk on FX Government bonds rate
rhofx4=0.077770133                 #Mark-down on FX Government bonds rate

#Rest of the World - Parameters-----------------------------------------------------------------------------------------------------------------------------------------------

sigmaRem=0.0002356039              #Share of World GDP distributed as remittances

varsigmafdi1=0.2                   #0.3 ->0.6 | #Speed of the shock on the FDI entering to the economy
varsigmafdi2=1.3                     #4 -> 1     |   #Initial period of the shock on the FDI entering to the economy 
varsigmafdi3=0.3665807     #Initial FDI entering to the economy as a share of NFCs investment 
varsigmafdi4=0.68                  #0.68 -> 0.60 |   #New FDI entering to the economy as a share of NFCs investment

zetaff=0.7762527                   #Share of FDI towards NFCs
shrGreenField=0.4840946            #Share of NFCs FDI that is greenfield

v1=0.008543549                     #Linear term in country risk premium
v2=2                               #Elasticity of country risk to the imports to official FX reserves ratio

zetabg=0.016                       #Share of the trade balance determining the purchase of domestic public bonds by the Rest of the World

betaen=2.5                         #Speed of convergence of nominal exchange rate

betaiwst=1                         #Speed of convergence of external short-term interest rate

alphaw=0.02                        #Foreign productivity growth rate
alphagw=0.03                       #World real GDP growth rate
alphapw=0.03                       #Foreign inflation rate
alphapO=0.03                       #Inflation rate of oil and coal prices


#Other - Parameters-----------------------------------------------------------------------------------------------------------------------------------------------

nuf = -0.0003904036                #Share of domestic production distributed as other flows to NFCs
nub = 0.007269213                  #Share of domestic production distributed as other flows to FCs
nug = 0.003295476                  #Share of domestic production distributed as other flows to the Government
nuh = 0.01509832                   #Share of domestic production distributed as other flows to Households
nuw = 0.00492403                   #Share of domestic production distributed as other flows to the Rest of the World


#Transition Scenario - Parameters-----------------------------------------------------------------------------------------------------------------------------------------------

shrDon = 0                                 #Share of NDC investment funded by international grants
shrGrTax = 0                               #Share of NDC investment funded by Green taxes
shrGrIC  = 0.5                             #Share of Green taxes levied on intermediate goods
shrGrL= 0                                  #Share of public investment funded with green loans
shrGrLFx = 0.5                             #Share of FX green loans
shrGrBw = 0                                #Share of domestic green bonds purchased by the RoW
md_lgtr = 0                              #Greenium on FX green loans
md_bgtr = 0                              #Greenium on domestic green bonds
K_0=1560.2222                              #NFCs' capital stock in 2019
lambdatr0 = 0                            #Speed of the NDC investment path
lambdatr1 = 15                            #Initial period of the NDC investment path
lambdatr2 = 0                             #Target NDC investment as a share of NFC's capital stock in 2019
taumtr=0.064085801                      #Import tax rate on transition goods (*Same as taum*)
tauothktr=0.024434906                   #Other tax rate on transition goods  (*Same as tauothk*)
betasigmamktr = 1.77651522362           #Speed of convergence of propensity to import transition goods (*Same as betasigmamk*)
reducXrO=-0.03                          #Reduction rate of fossil fuel exports
lambdatr0_adj = 0                            #Speed of the NDC investment path
lambdatr1_adj = 0                            #Initial period of the NDC investment path
lambdatr2_adj = 0                             #Target NDC investment as a share of NFC's capital stock in 2019

alpha_tr=0
beta_tr=0
gamma_tr=0
delta_tr=0

triggerFX = 0.12
speedFX = 80
shareFX = 0.5

speedTaum = 80
triggerTaum = 0.1381035
maxTaum = 6

speedWg = 80
triggerWG = 0.12
maxWg = 0

speedCg = 80
triggerCg = 0.12
maxCg = 0

triggerSTg = 0.12
speedSTg = 80
maxSTg = 0


#Alternative Scenario - Parameters----------------------------------------------------------------------------------------------------------------------------------------------
tauCBAM0=0                              #Linear term in CBAM tax rate
tauCBAM1=1                              #Sensitivity of CBAM tax rate to the ratio of "Green" to conventional NFCs' capital stock
sigmamktr0 = 0.58133                    #Autonomous propensity to import consumption goods
sigmaaktr = 0.00044                     #Linear term in productivity effect on propensity to import transition goods (*Double from sigmamk*)
epsilon2ktr = 0.241                     #Productivity elasticity for propensity to import transition goods  (*Same as epsilon2k*)
atr0=0                                  #Linear term in labour intensity of transition investments
atr1=1                                  #Sensitivity of labour intensity of transition investment to the ratio of "Green" to conventional NFCs' capital stock
scenInv=0                               #Government investment as a share of GDP to support exports diversification in the fossil fuels scenario



#Debt Swap Parameters
mdds=0
decds=0
betadebtSwapFXLgFX=4
betadebtSwapFXBgFX=4
betadebtSwapFXLgFXtr=4
dsactive=0
#Time---------------------------------------------------------------------------------------------------------------------------------------------------------------------------

##time

begin = 2019
end = 2050
by = 0.1
