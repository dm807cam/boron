###############################################
# Calculate the boron system from seawater pH #
# Dennis Mayk                                 # 
###############################################

pH2d11B <- function(date = "date", pH="pH", temp="temp", sal="sal", wd=0){
  
            if(length(pH) != length(temp) |
               length(pH) != length(sal)){
              stop('Vectors have different lengths.')
            }
  
            if(missing(date)) {
              pH2d11B_df <- data.frame(na.omit(cbind(pH, temp, sal, wd)))
            } else {
              pH2d11B_df <- data.frame(na.omit(cbind(date, pH, temp, sal, wd)))
            }
  
            pH <- pH2d11B_df$pH
            temp <- pH2d11B_df$temp
            sal <- pH2d11B_df$sal
            wd <- pH2d11B_df$wd
            
            # Constants
            Tk <- temp+273.15     # (K)
            Btot <- 432.6*sal/35  # (umol/kg), from Lee 2010
            P <- wd/10            # (bar = depth/10)
            H <- 10^(-1*pH)

            R <- 83.14472         # (cm^3.bar.mol^-1.K^-1)
            aB <- 1.0272
            eB <- (aB-1)*10^3
            d11Bsw <- 39.61       # (permil)
            R951 <- 4.04367
            RBsw <- (d11Bsw/1000+1)*R951

            # Kb calcuation after Dickson 1990
            tmp1 <- (-8966.9-2890.53*sqrt(sal)-77.942*sal+1.728*sal^(3/2)-0.0996*sal^2)
            tmp2 <- 148.0248+137.1942*sqrt(sal)+1.62142*sal
            tmp3 <- (-24.4344-25.085*sqrt(sal)-0.2474*sal)*log(temp+273.15)
            lnKb <-  tmp1/(temp+273.15) + tmp2 + tmp3 + 0.053105*sqrt(sal)*(temp+273.15)

            # T25,S35,Z0	8.597
            # T25,S35,Z3000	8.458

            # Kb pressure correction
            # formulation of Zeebe & Wolf-Gladrow (2001).  
            # The coefficients are often misquoted.  
            # The correct set are from Millero 1979 
            # (recalulated from Culberson and Pytlowicz 1968), 
            # as given by CO2SYS and Rae et al. (2011) and used here

            # lnkpok0(ipc) = -(deltav(ipc)./(R.*T)).*P + (0.5*deltak(ipc)./(R.*T)).*P.*P;
            # deltav(ipc)  =  a0(ipc) + a1(ipc).*TC + a2(ipc).*TC.*TC;

            a0 <- -29.48
            a1 <- 0.1622
            a2 <- -0.002608

            # deltak(ipc)  = (b0(ipc) + b1(ipc).*TC + b2(ipc).*TC.*TC); 

            b0 <- -0.00284
            b1 <- 0
            b2 <- 0
            Pfactor <- -(a0+a1*temp+a2*temp^2)/(R*Tk)*P+0.5*(b0+b1*temp+b2*temp^2)/(R*Tk)*P^2
            Kb_Pcor <- exp(lnKb+Pfactor)

            # Calculate 

            pH2d11B_df$`T (K)` <- Tk
            pH2d11B_df$`[B]tot (umol/kg), from Lee 2010` <- 432.6*sal/35
            pH2d11B_df$pKb_Pcor <- -log10(Kb_Pcor)   
            pH2d11B_df$`[B(OH)3]` <- pH2d11B_df$`[B]tot (umol/kg), from Lee 2010`/(1+Kb_Pcor/H)
            pH2d11B_df$`[B(OH)4-]` <- pH2d11B_df$`[B]tot (umol/kg), from Lee 2010`/(1+H/Kb_Pcor)
            pH2d11B_df$`R [B(OH)4-]` <- ((H^2*RBsw^2 + 
                                2*H^2*RBsw*aB + 
                                H^2*aB^2 + 
                                2*H*Kb_Pcor*RBsw^2*aB - 
                                2*H*Kb_Pcor*RBsw*aB^2 + 
                                8*H*Kb_Pcor*RBsw*aB - 
                                2*H*Kb_Pcor*RBsw + 
                                2*H*Kb_Pcor*aB + 
                                Kb_Pcor^2*RBsw^2*aB^2 + 
                                2*Kb_Pcor^2*RBsw*aB + 
                                Kb_Pcor^2)^(1/2) - 
                                H*aB - Kb_Pcor + 
                                H*RBsw + 
                                Kb_Pcor*RBsw*aB)/(2*aB*(H + Kb_Pcor))

            pH2d11B_df$`R [B(OH)3]` <- pH2d11B_df$`R [B(OH)4-]`*aB
            pH2d11B_df$`d11B (B(OH)4-)` <- (pH2d11B_df$`R [B(OH)4-]`/R951-1)*1000
            pH2d11B_df$`d11B (B(OH)3)` <- (pH2d11B_df$`R [B(OH)3]`/R951-1)*1000  
  
  return(pH2d11B_df)
            }

# Example 
pH2d11B(pH=8, temp=25, sal=35, wd=18)

