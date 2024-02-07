#' Timestamp generator
#'
#' Creating quarterly timestamps
#'
#' @param start_year Year of the beginning of the simulation. 
#' @param N_years Duration of the simulation in years.
#'
#' @return A vector of quarterly timestamps.
#' @examples
#'
#' @export
yyyy_mm_dd = function(start_year, N_years){
  grid = expand.grid(start_year:(start_year+N_years-1), c(1, 4, 7, 10))
  grid = grid[order(grid[,1], grid[,2]),]
  grid[,2] = as.character(grid[,2])
  grid[,2] = sapply(grid[,2], 
                    function(x){
                      ifelse(nchar(x)==1, 
                             paste0("0", x, collapse=''), 
                             x)
                    }
  )
  grid = rbind(grid, c(start_year+N_years, 1))
  times = as.Date(apply(grid, 1, 
                        function(x){
                          paste0(x[1], '-', x[2], '-01')
                        }
  )
  )
  times
}




#' Transition matrix
#'
#' Calculating the matrix of transition probabilities between states.
#'
#' @param p_hiv_infect The (marginal) probability of a random HIV negative individual to become HIV positive within a year. 
#' @param p_test The probability of HIV positive individuals getting tested (and diagnosed correctly).
#' @param p_ART The probability of HIV positive individuals deciding to go on ART.
#' @param p_suppress The probability of suppressed viral load among people on ART.
#' @param p_ART_adhere The probability of an individual on ART adhering to the treatment.
#' @param p_switch_therapy The probability of an ART patient with unsuppressed viral load choosing to switch to an alternative therapy.
#' @param p_ART_remorse The probability of individuals initially reluctant to undergo ART changing their mind.
#' @param p_death_gp Mortality rate of people who are HIV negative.
#' @param p_death_hiv Mortality rate of people who are HIV positive and are not on ART.
#' @param p_death_ART Mortality rate of people who are on ART with suppressed viral load.
#' @param p_death_ART_fail Mortality rate of people who are on ART with unsuppressed viral load
#'
#' @return A 5X6 stochastic matrix of annual transition probabilities.
#' @examples
#'
#' @export
trans_matrix = function(p_death_gp,
                        p_HIV_infect,
                        p_death_undiagnosed,
                        p_test,
                        p_death_pre_ART,
                        p_start_ART,
                        p_death_ART_suppressed,
                        p_ART_adhere,
                        p_ART_fail,
                        p_death_ART_nonsuppressed,
                        p_switch_ART,
                        p_alternative_ART_fail){
  P = matrix(0, nrow=6, ncol=6)
  P[1,6] = p_death_gp
  P[1,1:2] = (1 - p_death_gp)*c(1 - p_HIV_infect, p_HIV_infect)
  
  P[2,6] = p_death_undiagnosed
  P[2,2:3] = (1 - p_death_undiagnosed)*c(1 - p_test, p_test)
  
  P[3,6] = p_death_pre_ART
  P[3,3:4] = (1 - p_death_pre_ART)*c(1 - p_start_ART, p_start_ART)
  
  P[4,6] = p_death_ART_suppressed
  P[4,4:5] = (1 - p_death_ART_suppressed)*c(p_ART_adhere*(1 - p_ART_fail),
                                            1 - p_ART_adhere*(1 - p_ART_fail))
  
  P[5,6] = p_death_ART_nonsuppressed
  P[5,4:5] = (1 - p_death_ART_nonsuppressed)*c(p_switch_ART*(1 - p_alternative_ART_fail),
                                               1 - p_switch_ART*(1 - p_alternative_ART_fail))
  
  P[6,6] = 1
  
  P
}





#' HIV transmission probability
#'
#' The probability of contracting HIV within a year.
#'
#' @param p_partner_hiv The probability of a sexual partner being HIV positive. 
#' @param p_partner_ART_suppress The probability of a sexual partner being on ART with suppressed viral load.
#' @param p_condom The probability of using condom.
#' @param condom_effic Condom efficacy.
#' @param p_vaginal The probability of vaginal sex.
#' @param p_anal The probability of anal sex.
#' @param N_partners The average number of sexual partners per year.
#' @param acts_per_partner The average number of sexual acts per partner.
#' @param p_infected_vaginal_hiv The probability of getting infected through (unprotected) vaginal sex with HIV positive partners.
#' @param p_infected_anal_hiv The probability of getting infected through (unprotected) anal sex with HIV positive partners.
#' @param p_infected_vaginal_ART The probability of getting infected through (unprotected) vaginal sex with partners on ART with suppressed viral load.
#' @param p_infected_anal_ART The probability of getting infected through (unprotected) anal sex with partners on ART with suppressed viral load.
#' @param p_PrEP The probability of being on pre-exposure prophylaxis.
#' @param PrEP_effic The efficacy of being on pre-exposure prophylaxis.
#'
#' @return A real valued number between 0 and 1.
#' @examples
#'
#' @export
prob_infected = function(p_partner_hiv, p_partner_ART_suppress,
                         p_condom, condom_effic, p_vaginal, 
                         p_anal, N_partners, acts_per_partner,
                         p_infected_vaginal_hiv, p_infected_anal_hiv,
                         p_infected_vaginal_ART, p_infected_anal_ART, 
                         p_PrEP, PrEP_effic){
  coeff = (1 - condom_effic*p_condom)*(1 - PrEP_effic*p_PrEP)
  p_infected_in_act_hiv = coeff*(p_anal*p_infected_anal_hiv + 
                                   p_vaginal*p_infected_vaginal_hiv)
  p_infected_in_act_ART = coeff*(p_anal*p_infected_anal_ART + 
                                   p_vaginal*p_infected_vaginal_ART)
  a = (1-p_infected_in_act_hiv)^acts_per_partner*p_partner_hiv
  b = (1-p_infected_in_act_ART)^acts_per_partner*p_partner_ART_suppress
  p_partner_hiv_negative = 1 - p_partner_hiv - p_partner_ART_suppress
  
  1 - (a + b + p_partner_hiv_negative)^N_partners
}





#' HIV incidences
#'
#' Calculating incidences for the different populations groups.
#'
#' @param p_condom The probability of using condom among FSW.
#' @param p_condom_GP The probability of using condom among the general population.
#' @param condom_effic Condom efficacy.
#' @param p_vaginal The probability of vaginal sex.
#' @param p_anal The probability of anal sex.
#' @param N_partners_SW The average number of male partners each female sex worker meets per year.
#' @param acts_per_partner_SW The average number of sexual acts per male partner for a sex worker.
#' @param N_partners_GP The average number of sexual partners per year among the general public.
#' @param acts_per_partner_GP The average number of sexual acts per partner among the general public.
#' @param N_partners_MP The average number of different sex worker male partners meet yearly.
#' @param acts_per_partner_MP The average number of sexual acts per sex worker for a male partner.
#' @param p_infected_vaginal_hiv_SW The probability of getting infected through (unprotected) receptive vaginal sex with HIV positive partners.
#' @param p_infected_anal_hiv_SW The probability of getting infected through (unprotected) receptive anal sex with HIV positive partners.
#' @param p_infected_vaginal_ART_SW The probability of getting infected through (unprotected) receptive vaginal sex with partners on ART with suppressed viral load.
#' @param p_infected_anal_ART_SW The probability of getting infected through (unprotected) receptive anal sex with partners on ART with suppressed viral load.
#' @param p_infected_vaginal_hiv_MP The probability of getting infected through (unprotected) insertive vaginal sex with HIV positive partners.
#' @param p_infected_anal_hiv_MP The probability of getting infected through (unprotected) insertive anal sex with HIV positive partners.
#' @param p_infected_vaginal_ART_MP The probability of getting infected through (unprotected) insertive vaginal sex with partners on ART with suppressed viral load.
#' @param p_infected_anal_ART_MP The probability of getting infected through (unprotected) insertive anal sex with partners on ART with suppressed viral load.
#' @param MP_counts A vector of length 5 containing the number of male partners per (living) state.
#' @param SW_counts A vector of length 5 containing the number of female sex workers per (living) state.
#' @param GP_counts A vector of length 5 containing the number of pople from the general public per (living) state.
#' @param p_PrEP The probability of being on pre-exposure prophylaxis.
#' @param PrEP_effic The efficacy of being on pre-exposure prophylaxis.
#' @param MSM_counts A vector of length 5 containing the MSM count per (living) state.
#' @param N_partners_MSM The average number of sexual partners per year among MSM.
#' @param p_infected_anal_hiv_MSM The probability of getting infected through (unprotected) anal sex with HIV positive partners.
#' @param p_infected_anal_ART_MSM The probability of getting infected through (unprotected) anal sex with partners on ART.
#' @param acts_per_partner_MSM The average number of sexual acts per partner for MSM.
#' @param p_PrEP_MSM The probability of being on pre-exposure prophylaxis for MSM.
#' @param PrEP_effic_MSM The efficacy of being on pre-exposure prophylaxis for MSM.
#' @param p_condom_MSM The probability of using condom among MSM.
#'
#' @return A list containing the following items -
#' \item{prob_infected_SW}{Incidence for the female sex workers population.}
#' \item{prob_infected_MP}{Incidence for the male partners population.}
#' \item{prob_infected_GP}{Incidence for the general public.}
#' \item{prob_infected_MSM}{Incidence for MSM.}
#' @examples
#'
#' @export
Negative_to_Positive_prob = function(p_condom, p_condom_GP, 
                                     condom_effic, 
                                     p_vaginal, p_anal, 
                                     N_partners_SW, 
                                     acts_per_partner_SW,
                                     N_partners_GP,
                                     acts_per_partner_GP,
                                     N_partners_MP,
                                     acts_per_partner_MP,
                                     p_infected_vaginal_hiv_SW, 
                                     p_infected_anal_hiv_SW,
                                     p_infected_vaginal_ART_SW, 
                                     p_infected_anal_ART_SW,
                                     p_infected_vaginal_hiv_MP, 
                                     p_infected_anal_hiv_MP,
                                     p_infected_vaginal_ART_MP, 
                                     p_infected_anal_ART_MP,
                                     MP_counts, SW_counts, GP_counts,
                                     p_PrEP, PrEP_effic,
                                     MSM_counts, N_partners_MSM,
                                     p_infected_anal_hiv_MSM,
                                     p_infected_anal_ART_MSM,
                                     acts_per_partner_MSM,
                                     p_PrEP_MSM, PrEP_effic_MSM,
                                     p_condom_MSM){
  
  HIV_ART_props = function(counts){
    s = sum(counts)
    hiv_prop = sum(counts[c(2,3,5)])/s
    ART_prop = sum(counts[c(4)])/s
    
    list(hiv_prop=hiv_prop, ART_prop=ART_prop)
  }
  
  
  p_risk_MP = HIV_ART_props(MP_counts)
  prop_MP_hiv = p_risk_MP$hiv_prop
  prop_MP_ART = p_risk_MP$ART_prop
  
  prob_infected_SW = prob_infected(prop_MP_hiv, prop_MP_ART,
                                   p_condom, condom_effic, 
                                   p_vaginal, p_anal, 
                                   N_partners_SW, 
                                   acts_per_partner_SW,
                                   p_infected_vaginal_hiv_SW, 
                                   p_infected_anal_hiv_SW,
                                   p_infected_vaginal_ART_SW, 
                                   p_infected_anal_ART_SW,
                                   p_PrEP, PrEP_effic)
  
  p_risk_SW = HIV_ART_props(SW_counts)
  prop_SW_hiv = p_risk_SW$hiv_prop
  prop_SW_ART = p_risk_SW$ART_prop
  
  prob_infected_MP = prob_infected(prop_SW_hiv, prop_SW_ART,
                                   p_condom, condom_effic, 
                                   p_vaginal, p_anal, 
                                   N_partners_MP, 
                                   acts_per_partner_MP,
                                   p_infected_vaginal_hiv_MP, 
                                   p_infected_anal_hiv_MP,
                                   p_infected_vaginal_ART_MP, 
                                   p_infected_anal_ART_MP,
                                   p_PrEP=0, PrEP_effic)
  
  p_risk_GP = HIV_ART_props(GP_counts)
  prop_GP_hiv = p_risk_GP$hiv_prop
  prop_GP_ART = p_risk_GP$ART_prop
  
  prop_MP = sum(MP_counts)/(sum(MP_counts)+sum(GP_counts)*(1-p_female)+sum(MSM_counts))
  
  prob_infected_GP_MP = prob_infected(prop_MP_hiv, prop_MP_ART,
                                      p_condom_GP, condom_effic, 
                                      p_vaginal, p_anal, 
                                      N_partners_GP*prop_MP, 
                                      acts_per_partner_GP,
                                      p_infected_vaginal_hiv_MP, 
                                      p_infected_anal_hiv_MP,
                                      p_infected_vaginal_ART_MP, 
                                      p_infected_anal_ART_MP,
                                      p_PrEP=0, PrEP_effic)
  
  prob_infected_GP_GP_F = prob_infected(prop_GP_hiv, prop_GP_ART,
                                        p_condom_GP, condom_effic, 
                                        p_vaginal, p_anal, 
                                        N_partners_GP*(1-prop_MP), 
                                        acts_per_partner_GP,
                                        p_infected_vaginal_hiv_SW, 
                                        p_infected_anal_hiv_SW,
                                        p_infected_vaginal_ART_SW, 
                                        p_infected_anal_ART_SW,
                                        p_PrEP=0, PrEP_effic)
  
  prob_infected_GP_GP_M = prob_infected(prop_GP_hiv, prop_GP_ART,
                                        p_condom_GP, condom_effic, 
                                        p_vaginal, p_anal, 
                                        N_partners_GP, 
                                        acts_per_partner_GP,
                                        p_infected_vaginal_hiv_MP, 
                                        p_infected_anal_hiv_MP,
                                        p_infected_vaginal_ART_MP, 
                                        p_infected_anal_ART_MP,
                                        p_PrEP=0, PrEP_effic)
  
  
  
  p_risk_MSM = HIV_ART_props(MSM_counts)
  prop_MSM_hiv = p_risk_MSM$hiv_prop
  prop_MSM_ART = p_risk_MSM$ART_prop
  
  prop_MSM = sum(MSM_counts)/(sum(MP_counts)+sum(GP_counts)*(1-p_female)+sum(MSM_counts))
  
  prob_infected_GP_MSM = prob_infected(prop_MSM_hiv, prop_MSM_ART,
                                       p_condom_GP, condom_effic, 
                                       p_vaginal, p_anal, 
                                       N_partners_GP*prop_MSM, 
                                       acts_per_partner_GP,
                                       p_infected_vaginal_hiv_MP, 
                                       p_infected_anal_hiv_MSM,
                                       p_infected_vaginal_ART_MP, 
                                       p_infected_anal_ART_MSM,
                                       p_PrEP_MSM, PrEP_effic_MSM)
  
  
  
  prob_infected_MSM = prob_infected(prop_MSM_hiv, prop_MSM_ART,
                                    p_condom_MSM, condom_effic, 
                                    p_vaginal=0, p_anal=1, 
                                    N_partners_MSM, 
                                    acts_per_partner_MSM,
                                    p_infected_vaginal_hiv_MP, 
                                    p_infected_anal_hiv_MSM,
                                    p_infected_vaginal_ART_MP, 
                                    p_infected_anal_ART_MSM,
                                    p_PrEP_MSM, PrEP_effic_MSM)
  
  prob_infected_GP = prob_infected_GP_MP*prop_MP + 
    prob_infected_GP_GP_F*(1-prop_MP-prop_MSM)*p_female +
    prob_infected_GP_GP_M*(1-prop_MP-prop_MSM)*(1-p_female) +
    prob_infected_GP_MSM*prop_MSM
  
  
  list(prob_infected_SW=prob_infected_SW, 
       prob_infected_MP=prob_infected_MP,
       prob_infected_GP=prob_infected_GP,
       prob_infected_MSM=prob_infected_MSM)
}





#' Quarterly demographic update
#'
#' Updating the counts of the different groups across the different states.
#'
#' @param p_condom The probability of using condom.
#' @param p_condom_GP The probability of using condom among the general population.
#' @param condom_effic Condom efficacy.
#' @param p_vaginal The probability of vaginal sex.
#' @param p_anal The probability of anal sex.
#' @param N_partners_GP The average number of sexual partners per year among the general public.
#' @param acts_per_partner_GP The average number of sexual acts per partner among the general public.
#' @param N_partners_MP The average number of different sex worker male partners meet yearly.
#' @param acts_per_partner_MP The average number of sexual acts per sex worker for a male partner.
#' @param p_infected_vaginal_hiv_SW The probability of getting infected through (unprotected) receptive vaginal sex with HIV positive partners.
#' @param p_infected_anal_hiv_SW The probability of getting infected through (unprotected) receptive anal sex with HIV positive partners.
#' @param p_infected_vaginal_ART_SW The probability of getting infected through (unprotected) receptive vaginal sex with partners on ART with suppressed viral load.
#' @param p_infected_anal_ART_SW The probability of getting infected through (unprotected) receptive anal sex with partners on ART with suppressed viral load.
#' @param p_infected_vaginal_hiv_MP The probability of getting infected through (unprotected) insertive vaginal sex with HIV positive partners.
#' @param p_infected_anal_hiv_MP The probability of getting infected through (unprotected) insertive anal sex with HIV positive partners.
#' @param p_infected_vaginal_ART_MP The probability of getting infected through (unprotected) insertive vaginal sex with partners on ART with suppressed viral load.
#' @param p_infected_anal_ART_MP The probability of getting infected through (unprotected) insertive anal sex with partners on ART with suppressed viral load.
#' @param MP_counts A vector of length 5 containing the number of male partners per (living) state.
#' @param SW_counts A vector of length 5 containing the number of female sex workers per (living) state.
#' @param GP_counts A vector of length 5 containing the number of pople from the general public per (living) state.
#' @param p_test_SW The probability of HIV positive sex workers getting tested (and diagnosed correctly).
#' @param p_test_GP The probability of HIV positive individuals in the general population getting tested (and diagnosed correctly).
#' @param p_ART The probability of HIV positive individuals deciding to go on ART.
#' @param p_suppress The probability of suppressed viral load among people on ART.
#' @param p_ART_adhere The probability of an individual on ART adhering to the treatment.
#' @param p_switch_therapy The probability of an ART patient with unsuppressed viral load choosing to switch to an alternative therapy.
#' @param p_ART_remorse The probability of individuals initially reluctant to undergo ART changing their mind.
#' @param p_death_gp Mortality rate of people who are HIV negative.
#' @param p_death_hiv Mortality rate of people who are HIV positive and are not on ART.
#' @param p_death_ART Mortality rate of people who are on ART with suppressed viral load.
#' @param p_death_ART_fail Mortality rate of people who are on ART with unsuppressed viral load
#' @param p_PrEP The probability of being on pre-exposure prophylaxis.
#' @param PrEP_effic The efficacy of being on pre-exposure prophylaxis.
#' @param SW_retire_rate Proportion of SWs leaving their population (to join the GP) quarterly.
#' @param SW_recruit_rate Proportion of females in the GP becoming SWs quarterly.
#' @param MP_retire_rate Proportion of MPs leaving their population (to join the GP) quarterly.
#' @param MP_recruit_rate Proportion of males in the GP becoming MPs quarterly.
#' @param growth_rate Quarterly growth rate in the Rwandan population.
#' @param MSM_counts A vector of length 5 containing the MSM count per (living) state.
#' @param N_partners_MSM The average number of sexual partners per year among MSM.
#' @param p_infected_anal_hiv_MSM The probability of getting infected through (unprotected) anal sex with HIV positive partners.
#' @param p_infected_anal_ART_MSM The probability of getting infected through (unprotected) anal sex with partners on ART.
#' @param acts_per_partner_MSM The average number of sexual acts per partner for MSM.
#' @param p_PrEP_MSM The probability of being on pre-exposure prophylaxis for MSM.
#' @param PrEP_effic_MSM The efficacy of being on pre-exposure prophylaxis for MSM.
#' @param p_condom_MSM The probability of using condom among MSM.
#' @param p_test_MSM The probability of HIV positive MSM getting tested (and diagnosed correctly).
#' @param MSM_retire_rate Proportion of MSM leaving their population (to join the GP) quarterly.
#' @param MSM_recruit_rate Proportion of males in the GP becoming MSM quarterly.
#'
#' @return A list containing the following items -
#' \item{SW_counts}{A vector of length 5 containing the updated number of female sex workers per (living) state.}
#' \item{MP_counts}{A vector of length 5 containing the updated number of male partners per (living) state.}
#' \item{GP_counts}{A vector of length 5 containing the updated number of people among the general public per (living) state.}
#' \item{MSM_counts}{A vector of length 5 containing the updated number of people among the MSM per (living) state.}
#' \item{incidences}{A vector of length 4 containing the incidences of the four groups (SW, MP, MSM and GP).}
#' \item{overall_incidence}{The overall incidence of HIV for the entire population.}
#' @examples
#'
#' @export
Quarterly_Cycle = function(p_condom, p_condom_GP,
                           condom_effic, 
                           p_vaginal, p_anal, 
                           N_partners_GP,
                           acts_per_partner_GP,
                           N_partners_MP,
                           acts_per_partner_MP,
                           p_infected_vaginal_hiv_SW, 
                           p_infected_anal_hiv_SW,
                           p_infected_vaginal_ART_SW, 
                           p_infected_anal_ART_SW,
                           p_infected_vaginal_hiv_MP, 
                           p_infected_anal_hiv_MP,
                           p_infected_vaginal_ART_MP, 
                           p_infected_anal_ART_MP,
                           MP_counts, SW_counts, GP_counts,
                           p_test_SW, p_test_GP,
                           p_death_pre_ART,
                           p_start_ART,
                           p_death_ART_suppressed,
                           p_ART_adhere,
                           p_ART_fail,
                           p_death_ART_nonsuppressed,
                           p_switch_ART,
                           p_alternative_ART_fail,
                           p_PrEP, PrEP_effic,
                           SW_retire_rate, SW_recruit_rate,
                           MP_retire_rate, MP_recruit_rate,
                           growth_rate,
                           MSM_counts, N_partners_MSM,
                           p_infected_anal_hiv_MSM,
                           p_infected_anal_ART_MSM,
                           acts_per_partner_MSM,
                           p_PrEP_MSM, PrEP_effic_MSM,
                           p_condom_MSM, p_test_MSM,
                           MSM_retire_rate, MSM_recruit_rate){
  acts_per_partner_SW = acts_per_partner_MP
  N_partners_SW = N_partners_MP*acts_per_partner_MP*sum(MP_counts)
  N_partners_SW = N_partners_SW/sum(SW_counts)/acts_per_partner_SW
  
  
  infect_probs = Negative_to_Positive_prob(p_condom, p_condom_GP, 
                                           condom_effic, 
                                           p_vaginal, p_anal, 
                                           N_partners_SW, 
                                           acts_per_partner_SW,
                                           N_partners_GP,
                                           acts_per_partner_GP,
                                           N_partners_MP,
                                           acts_per_partner_MP,
                                           p_infected_vaginal_hiv_SW, 
                                           p_infected_anal_hiv_SW,
                                           p_infected_vaginal_ART_SW, 
                                           p_infected_anal_ART_SW,
                                           p_infected_vaginal_hiv_MP, 
                                           p_infected_anal_hiv_MP,
                                           p_infected_vaginal_ART_MP, 
                                           p_infected_anal_ART_MP,
                                           MP_counts, SW_counts, GP_counts,
                                           p_PrEP, PrEP_effic,
                                           MSM_counts, N_partners_MSM,
                                           p_infected_anal_hiv_MSM,
                                           p_infected_anal_ART_MSM,
                                           acts_per_partner_MSM,
                                           p_PrEP_MSM, PrEP_effic_MSM,
                                           p_condom_MSM)
  
  p_hiv_infected_SW = infect_probs$prob_infected_SW
  p_hiv_infected_MP = infect_probs$prob_infected_MP
  p_hiv_infected_GP = infect_probs$prob_infected_GP
  p_hiv_infected_MSM = infect_probs$prob_infected_MSM
  
  P_trans_SW = trans_matrix(p_death_gp,
                            p_hiv_infected_SW,
                            p_death_undiagnosed,
                            p_test_SW,
                            p_death_pre_ART,
                            p_start_ART,
                            p_death_ART_suppressed,
                            p_ART_adhere,
                            p_ART_fail,
                            p_death_ART_nonsuppressed,
                            p_switch_ART,
                            p_alternative_ART_fail)
  
  P_trans_MP = trans_matrix(p_death_gp,
                            p_hiv_infected_MP,
                            p_death_undiagnosed,
                            p_test_GP,
                            p_death_pre_ART,
                            p_start_ART,
                            p_death_ART_suppressed,
                            p_ART_adhere,
                            p_ART_fail,
                            p_death_ART_nonsuppressed,
                            p_switch_ART,
                            p_alternative_ART_fail)
  
  P_trans_GP = trans_matrix(p_death_gp,
                            p_hiv_infected_GP,
                            p_death_undiagnosed,
                            p_test_GP,
                            p_death_pre_ART,
                            p_start_ART,
                            p_death_ART_suppressed,
                            p_ART_adhere,
                            p_ART_fail,
                            p_death_ART_nonsuppressed,
                            p_switch_ART,
                            p_alternative_ART_fail)
  
  P_trans_MSM = trans_matrix(p_death_gp,
                             p_hiv_infected_MSM,
                             p_death_undiagnosed,
                             p_test_MSM,
                             p_death_pre_ART,
                             p_start_ART,
                             p_death_ART_suppressed,
                             p_ART_adhere,
                             p_ART_fail,
                             p_death_ART_nonsuppressed,
                             p_switch_ART,
                             p_alternative_ART_fail)
  
  
  incidences = c(p_hiv_infected_GP, p_hiv_infected_SW, 
                 p_hiv_infected_MP, p_hiv_infected_MSM)
  counts = c(sum(GP_counts), sum(SW_counts), 
             sum(MP_counts), sum(MSM_counts))
  overall_incidence = incidences%*%counts/sum(counts)
  
  redistribute = function(counts, trans_probs){
    aux_func = function(i){
      rmultinom(1, counts[i], trans_probs[i,])
    }
    
    apply(sapply(1:length(counts), aux_func), 1, sum)
  }
  
  
  Old_population = sum(c(GP_counts, SW_counts, MP_counts, MSM_counts))
  SW_counts_new = redistribute(SW_counts, P_trans_SW)
  MP_counts_new = redistribute(MP_counts, P_trans_MP)
  GP_counts_new = redistribute(GP_counts, P_trans_GP)
  MSM_counts_new = redistribute(MSM_counts, P_trans_MSM)
  
  Deaths = sum(SW_counts_new[6], MP_counts_new[6], 
               GP_counts_new[6], MSM_counts_new[6])
  Births = round(Old_population*growth_rate + Deaths)
  SW_counts_new = SW_counts_new[1:5]
  MP_counts_new = MP_counts_new[1:5]
  MSM_counts_new = MSM_counts_new[1:5]
  GP_counts_new = GP_counts_new[1:5]
  GP_counts_new[1] = GP_counts_new[1] + Births
  SW_retirement = redistribute(round(SW_retire_rate*sum(SW_counts)), 
                               t(as.matrix(SW_counts/sum(SW_counts))))
  SW_recruitment = redistribute(round(SW_recruit_rate*sum(GP_counts)*sum(SW_counts)/Old_population), 
                                t(as.matrix(GP_counts/sum(GP_counts))))
  MP_retirement = redistribute(round(MP_retire_rate*sum(MP_counts)), 
                               t(as.matrix(MP_counts/sum(MP_counts))))
  MP_recruitment = redistribute(round(MP_recruit_rate*sum(GP_counts)*sum(MP_counts)/Old_population), 
                                t(as.matrix(GP_counts/sum(GP_counts))))
  MSM_retirement = redistribute(round(MSM_retire_rate*sum(MSM_counts)), 
                                t(as.matrix(GP_counts/sum(GP_counts))))
  MSM_recruitment = redistribute(round(MSM_recruit_rate*sum(GP_counts)*sum(MSM_counts)/Old_population), 
                                 t(as.matrix(GP_counts/sum(GP_counts))))
  
  GP_counts = GP_counts_new + MP_retirement - MP_recruitment + 
    SW_retirement - SW_recruitment + MSM_retirement - MSM_recruitment
  MP_counts = MP_counts_new + MP_recruitment - MP_retirement
  SW_counts = SW_counts_new + SW_recruitment - SW_retirement
  MSM_counts = MSM_counts_new + MSM_recruitment - MSM_retirement
  
  return(list(SW_counts=SW_counts, 
              MP_counts=MP_counts, 
              GP_counts=GP_counts,
              MSM_counts=MSM_counts,
              incidences=1 - (1 - incidences)^4,
              overall_incidence=1 - (1 - overall_incidence)^4)
  )
}






#' A single simulation
#'
#' Simulating once over a fixed time period.
#'
#' @param N_years Duration of the simulation in years.
#' @param p_condom_GP The probability of using condom among the general population.
#' @param p_condom_low The probability of using condom among SW at the beginning of the simulation.
#' @param p_condom_high The probability of using condom among SW at the end of the simulation.
#' @param condom_effic Condom efficacy.
#' @param p_vaginal The probability of vaginal sex.
#' @param p_anal The probability of anal sex.
#' @param N_partners_GP The average number of sexual partners per year among the general public.
#' @param acts_per_partner_GP The average number of sexual acts per partner among the general public.
#' @param N_partners_MP The average number of different sex worker male partners meet yearly.
#' @param acts_per_partner_MP The average number of sexual acts per sex worker for a male partner.
#' @param p_infected_vaginal_hiv_SW The probability of getting infected through (unprotected) receptive vaginal sex with HIV positive partners.
#' @param p_infected_anal_hiv_SW The probability of getting infected through (unprotected) receptive anal sex with HIV positive partners.
#' @param p_infected_vaginal_ART_SW The probability of getting infected through (unprotected) receptive vaginal sex with partners on ART with suppressed viral load.
#' @param p_infected_anal_ART_SW The probability of getting infected through (unprotected) receptive anal sex with partners on ART with suppressed viral load.
#' @param p_infected_vaginal_hiv_MP The probability of getting infected through (unprotected) insertive vaginal sex with HIV positive partners.
#' @param p_infected_anal_hiv_MP The probability of getting infected through (unprotected) insertive anal sex with HIV positive partners.
#' @param p_infected_vaginal_ART_MP The probability of getting infected through (unprotected) insertive vaginal sex with partners on ART with suppressed viral load.
#' @param p_infected_anal_ART_MP The probability of getting infected through (unprotected) insertive anal sex with partners on ART with suppressed viral load.
#' @param MP_counts A vector of length 5 containing the number of male partners per (living) state.
#' @param SW_counts A vector of length 5 containing the number of female sex workers per (living) state.
#' @param GP_counts A vector of length 5 containing the number of pople from the general public per (living) state.
#' @param p_test_SW The probability of HIV positive sex workers getting tested (and diagnosed correctly).
#' @param p_test_GP The probability of HIV positive individuals in the general population getting tested (and diagnosed correctly).
#' @param p_ART The probability of HIV positive individuals deciding to go on ART.
#' @param p_suppress The probability of suppressed viral load among people on ART.
#' @param p_ART_adhere_low The probability of an individual on ART adhering to the treatment at the beginning of the simulation.
#' @param p_ART_adhere_high The probability of an individual on ART adhering to the treatment at the end of the simulation.
#' @param p_switch_therapy The probability of an ART patient with unsuppressed viral load choosing to switch to an alternative therapy.
#' @param p_ART_remorse The probability of individuals initially reluctant to undergo ART changing their mind.
#' @param p_death_gp Mortality rate of people who are HIV negative.
#' @param p_death_hiv Mortality rate of people who are HIV positive and are not on ART.
#' @param p_death_ART Mortality rate of people who are on ART with suppressed viral load.
#' @param p_death_ART_fail Mortality rate of people who are on ART with unsuppressed viral load
#' @param p_PrEP_low The probability of a SW being on pre-exposure prophylaxis at the beginning of the simulation.
#' @param p_PrEP_up The probability of a SW being on pre-exposure prophylaxis at the end of the simulation.
#' @param PrEP_effic The efficacy of being on pre-exposure prophylaxis for SW.
#' @param time_to_PrEP Time (in years) from beginning of simulation to introduction of PrEP.
#' @param SW_retire_rate Proportion of SWs leaving their population (to join the GP) quarterly.
#' @param SW_recruit_rate Proportion of females in the GP becoming SWs quarterly.
#' @param MP_retire_rate Proportion of MPs leaving their population (to join the GP) quarterly.
#' @param MP_recruit_rate Proportion of males in the GP becoming MPs quarterly.
#' @param growth_rate Quarterly growth rate in the Rwandan population.
#' @param MSM_counts A vector of length 5 containing the MSM count per (living) state.
#' @param N_partners_MSM The average number of sexual partners per year among MSM.
#' @param p_infected_anal_hiv_MSM The probability of getting infected through (unprotected) anal sex with HIV positive partners.
#' @param p_infected_anal_ART_MSM The probability of getting infected through (unprotected) anal sex with partners on ART.
#' @param acts_per_partner_MSM The average number of sexual acts per partner for MSM.
#' @param p_PrEP_MSM_low The probability of being on pre-exposure prophylaxis for MSM at the beginning of the simulation.
#' @param p_PrEP_MSM_up The probability of being on pre-exposure prophylaxis for MSM at the end of the simulation.
#' @param PrEP_effic_MSM The efficacy of being on pre-exposure prophylaxis for MSM.
#' @param p_condom_MSM_low The probability of using condom among MSM at the beginning of the simulation.
#' @param p_condom_MSM_high The probability of using condom among MSM at the end of the simulation.
#'
#' @return A list containing the following items -
#' \item{SW_counts}{A vector of length 5 containing the updated number of female sex workers per (living) state.}
#' \item{MP_counts}{A vector of length 5 containing the updated number of male partners per (living) state.}
#' \item{GP_counts}{A vector of length 5 containing the updated number of people among the general public per (living) state.}
#' \item{MSM_counts}{A vector of length 5 containing the updated number of people among the MSM per (living) state.}
#' \item{incidences}{A vector of length 4 containing the incidences of the four groups (SW, MP, MSM and GP).}
#' \item{overall_incidence}{The overall incidence of HIV for the entire population.}
#' @examples
#'
#' @export
single_simulation = function(N_years, p_condom_GP,
                             p_condom_low, p_condom_high, 
                             condom_effic, p_vaginal, p_anal, 
                             N_partners_GP, acts_per_partner_GP, 
                             N_partners_MP, acts_per_partner_MP,
                             p_infected_vaginal_hiv_SW, 
                             p_infected_anal_hiv_SW,
                             p_infected_vaginal_ART_SW, 
                             p_infected_anal_ART_SW,
                             p_infected_vaginal_hiv_MP, 
                             p_infected_anal_hiv_MP,
                             p_infected_vaginal_ART_MP, 
                             p_infected_anal_ART_MP,
                             MP_counts, SW_counts, GP_counts,
                             p_test_SW, p_test_GP,
                             p_death_pre_ART,
                             p_start_ART_low,
                             p_start_ART_high,
                             p_death_ART_suppressed,
                             p_ART_adhere_low,
                             p_ART_adhere_high,
                             p_ART_fail,
                             p_death_ART_nonsuppressed,
                             p_switch_ART,
                             p_alternative_ART_fail,
                             p_PrEP_low, p_PrEP_up,
                             PrEP_effic, time_to_PrEP,
                             SW_retire_rate, SW_recruit_rate,
                             MP_retire_rate, MP_recruit_rate,
                             growth_rate,
                             MSM_counts, N_partners_MSM,
                             p_infected_anal_hiv_MSM,
                             p_infected_anal_ART_MSM,
                             acts_per_partner_MSM,
                             p_PrEP_MSM_low, p_PrEP_MSM_up,
                             PrEP_effic_MSM, p_condom_MSM_low, 
                             p_condom_MSM_high, p_test_MSM,
                             MSM_retire_rate, MSM_recruit_rate){
  GP_matrix = MP_matrix = SW_matrix = MSM_matrix = matrix(0, nrow=N_years*4+1, ncol=5)
  incidence_matrix = matrix(0, nrow=N_years*4+1, ncol=4)
  overall_incidence = rep(0, N_years*4+1)
  GP_matrix[1,] = GP_counts
  SW_matrix[1,] = SW_counts
  MP_matrix[1,] = MP_counts
  MSM_matrix[1,] = MSM_counts
  
  year_PrEP_intro = 1 + time_to_PrEP
  delta = (p_PrEP_up - p_PrEP_low)/(N_years-time_to_PrEP)/4
  delta_MSM = (p_PrEP_MSM_up - p_PrEP_MSM_low)/(N_years-time_to_PrEP)/4
  a_exp = p_start_ART_low*(p_start_ART_low/p_start_ART_high)^(time_to_PrEP/N_years)
  b_exp = log(p_start_ART_low/p_start_ART_high)/N_years
  b_log = N_years/(exp(p_start_ART_high-p_start_ART_low)-1)
  a_log = p_start_ART_low - log(b_log)
  delta_ART = (p_start_ART_high - p_start_ART_low)/(N_years*4+1)
  delta_condom = (p_condom_high - p_condom_low)/(N_years*4+1)
  delta_condom_MSM = (p_condom_MSM_high - p_condom_MSM_low)/(N_years*4+1)
  delta_ART_adhere = (p_ART_adhere_high - p_ART_adhere_low)/(N_years*4+1)
  
  for(i in 5:(4*N_years+6)){
    prev_PrEP = ifelse(i/4 <= year_PrEP_intro, 0, 
                       p_PrEP_low + delta*(i-4*time_to_PrEP-1))
    prev_PrEP_MSM = ifelse(i/4 <= year_PrEP_intro, 0, 
                           p_PrEP_MSM_low + delta_MSM*(i-4*time_to_PrEP-1))
    
    p_ART = ifelse(p_start_ART_low < p_start_ART_high,
                   a_log + log(b_log + (i-5)/4),
                   p_start_ART_low)
    p_condom = p_condom_low + delta_condom*(i-5)
    p_condom_MSM = p_condom_MSM_low + delta_condom_MSM*(i-5)
    p_ART_adhere = p_ART_adhere_low + delta_ART_adhere*(i-5)
    population = Quarterly_Cycle(p_condom, p_condom_GP, 
                                 condom_effic, p_vaginal, p_anal, 
                                 N_partners_GP,
                                 acts_per_partner_GP,
                                 N_partners_MP,
                                 acts_per_partner_MP,
                                 p_infected_vaginal_hiv_SW, 
                                 p_infected_anal_hiv_SW,
                                 p_infected_vaginal_ART_SW, 
                                 p_infected_anal_ART_SW,
                                 p_infected_vaginal_hiv_MP, 
                                 p_infected_anal_hiv_MP,
                                 p_infected_vaginal_ART_MP, 
                                 p_infected_anal_ART_MP,
                                 MP_counts, SW_counts, GP_counts,
                                 p_test_SW, p_test_GP,
                                 p_death_pre_ART,
                                 p_ART,
                                 p_death_ART_suppressed,
                                 p_ART_adhere,
                                 p_ART_fail,
                                 p_death_ART_nonsuppressed,
                                 p_switch_ART,
                                 p_alternative_ART_fail,
                                 prev_PrEP, PrEP_effic,
                                 SW_retire_rate, SW_recruit_rate,
                                 MP_retire_rate, MP_recruit_rate,
                                 growth_rate,
                                 MSM_counts, N_partners_MSM,
                                 p_infected_anal_hiv_MSM,
                                 p_infected_anal_ART_MSM,
                                 acts_per_partner_MSM,
                                 prev_PrEP_MSM, PrEP_effic_MSM,
                                 p_condom_MSM, p_test_MSM,
                                 MSM_retire_rate, MSM_recruit_rate)
    
    GP_counts = population$GP_counts
    SW_counts = population$SW_counts
    MP_counts = population$MP_counts
    MSM_counts = population$MSM_counts
    
    overall_incidence[i-5] = population$overall_incidence
    incidence_matrix[i-5,] = population$incidences
    
    GP_matrix[i-5,] = GP_counts
    SW_matrix[i-5,] = SW_counts
    MP_matrix[i-5,] = MP_counts
    MSM_matrix[i-5,] = MSM_counts
  }
  
  HIV_GP = apply(GP_matrix[,-1], 1, sum)
  HIV_SW = apply(SW_matrix[,-1], 1, sum)
  HIV_MP = apply(MP_matrix[,-1], 1, sum)
  HIV_MSM = apply(MSM_matrix[,-1], 1, sum)
  HIV = HIV_GP + HIV_SW + HIV_MP + HIV_MSM
  
  N_GP = apply(GP_matrix, 1, sum)
  N_SW = apply(SW_matrix, 1, sum)
  N_MP = apply(MP_matrix, 1, sum)
  N_MSM = apply(MSM_matrix, 1, sum)
  Total_pop = N_GP + N_SW + N_MP + N_MSM
  
  HIV_GP_prop = HIV_GP/N_GP*100
  HIV_SW_prop = HIV_SW/N_SW*100
  HIV_MP_prop = HIV_MP/N_MP*100
  HIV_MSM_prop = HIV_MSM/N_MSM*100
  HIV_prop = HIV/Total_pop*100
  
  ART_GP = apply(GP_matrix[,4:5], 1, sum)
  ART_SW = apply(SW_matrix[,4:5], 1, sum)
  ART_MP = apply(MP_matrix[,4:5], 1, sum)
  ART_MSM = apply(MSM_matrix[,4:5], 1, sum)
  ART = ART_GP + ART_SW + ART_MP + ART_MSM
  
  ART_GP_prop = ART_GP/HIV_GP*100
  ART_SW_prop = ART_SW/HIV_SW*100
  ART_MP_prop = ART_MP/HIV_MP*100
  ART_MSM_prop = ART_MSM/HIV_MSM*100
  ART_prop = ART/HIV*100
  
  return(list(GP=GP_matrix, SW=SW_matrix, 
              MP=MP_matrix, MSM=MSM_matrix,
              incidences=incidence_matrix, 
              overall_incidence=overall_incidence,
              HIV_counts=HIV, HIV_props=HIV_prop, 
              ART_counts=ART, ART_props=ART_prop,
              HIV_GP=HIV_GP, HIV_SW=HIV_SW, 
              HIV_MP=HIV_MP, HIV_MSM=HIV_MSM,
              ART_GP=ART_GP, ART_SW=ART_SW, 
              ART_MP=ART_MP, ART_MSM=ART_MSM,
              HIV_GP_prop=HIV_GP_prop, HIV_SW_prop=HIV_SW_prop, 
              HIV_MP_prop=HIV_MP_prop, HIV_MSM_prop=HIV_MSM_prop,
              ART_GP_prop=ART_GP_prop, ART_SW_prop=ART_SW_prop, 
              ART_MP_prop=ART_MP_prop, ART_MSM_prop=ART_MSM_prop))
}





simulation = function(N_simulations, N_years, p_condom_GP,
                      p_condom_low, p_condom_high, 
                      condom_effic, p_vaginal, 
                      p_anal, N_partners_GP,
                      acts_per_partner_GP, N_partners_MP,
                      acts_per_partner_MP,
                      p_infected_vaginal_hiv_SW, 
                      p_infected_anal_hiv_SW,
                      p_infected_vaginal_ART_SW, 
                      p_infected_anal_ART_SW,
                      p_infected_vaginal_hiv_MP, 
                      p_infected_anal_hiv_MP,
                      p_infected_vaginal_ART_MP, 
                      p_infected_anal_ART_MP,
                      MP_counts, SW_counts, GP_counts,
                      p_test_SW, p_test_GP,
                      p_death_pre_ART,
                      p_start_ART_low,
                      p_start_ART_high,
                      p_death_ART_suppressed,
                      p_ART_adhere_low,
                      p_ART_adhere_high,
                      p_ART_fail,
                      p_death_ART_nonsuppressed,
                      p_switch_ART,
                      p_alternative_ART_fail,
                      p_PrEP_low, p_PrEP_up,
                      PrEP_effic, time_to_PrEP,
                      SW_retire_rate, SW_recruit_rate,
                      MP_retire_rate, MP_recruit_rate,
                      growth_rate,
                      MSM_counts, N_partners_MSM,
                      p_infected_anal_hiv_MSM,
                      p_infected_anal_ART_MSM,
                      acts_per_partner_MSM,
                      p_PrEP_MSM_low, p_PrEP_MSM_up,
                      PrEP_effic_MSM, p_condom_MSM_low, 
                      p_condom_MSM_high, p_test_MSM,
                      MSM_retire_rate, MSM_recruit_rate){
  sim_list = list()
  for(i in 1:N_simulations){
    sim_list = c(sim_list, 
                 single_simulation(N_years, 
                                   p_condom_GP, p_condom_low, 
                                   p_condom_high, condom_effic, 
                                   p_vaginal, p_anal, N_partners_GP,
                                   acts_per_partner_GP, N_partners_MP,
                                   acts_per_partner_MP,
                                   p_infected_vaginal_hiv_SW, 
                                   p_infected_anal_hiv_SW,
                                   p_infected_vaginal_ART_SW, 
                                   p_infected_anal_ART_SW,
                                   p_infected_vaginal_hiv_MP, 
                                   p_infected_anal_hiv_MP,
                                   p_infected_vaginal_ART_MP, 
                                   p_infected_anal_ART_MP,
                                   MP_counts, SW_counts, GP_counts,
                                   p_test_SW, p_test_GP,
                                   p_death_pre_ART,
                                   p_start_ART_low,
                                   p_start_ART_high,
                                   p_death_ART_suppressed,
                                   p_ART_adhere_low,
                                   p_ART_adhere_high,
                                   p_ART_fail,
                                   p_death_ART_nonsuppressed,
                                   p_switch_ART,
                                   p_alternative_ART_fail,
                                   p_PrEP_low, p_PrEP_up,
                                   PrEP_effic, time_to_PrEP,
                                   SW_retire_rate, SW_recruit_rate,
                                   MP_retire_rate, MP_recruit_rate,
                                   growth_rate,
                                   MSM_counts, N_partners_MSM,
                                   p_infected_anal_hiv_MSM,
                                   p_infected_anal_ART_MSM,
                                   acts_per_partner_MSM,
                                   p_PrEP_MSM_low, p_PrEP_MSM_up,
                                   PrEP_effic_MSM, p_condom_MSM_low, 
                                   p_condom_MSM_high, p_test_MSM,
                                   MSM_retire_rate, MSM_recruit_rate))
  }
  sim_list
}




simulation_analysis = function(sim, alpha, start_year, N_simulations){
  GPs = sim[26*(0:(N_simulations-1))+1]
  SWs = sim[26*(0:(N_simulations-1))+2]
  MPs = sim[26*(0:(N_simulations-1))+3]
  MSMs = sim[26*(0:(N_simulations-1))+4]
  incidences = sim[26*(0:(N_simulations-1))+5]
  
  overall_incidences = 1000*do.call(rbind, sim[26*(0:(N_simulations-1))+6])
  HIV_counts = do.call(rbind, sim[26*(0:(N_simulations-1))+7])
  HIV_props = do.call(rbind, sim[26*(0:(N_simulations-1))+8])
  ART_counts = do.call(rbind, sim[26*(0:(N_simulations-1))+9])
  ART_props = do.call(rbind, sim[26*(0:(N_simulations-1))+10])
  HIV_GP = do.call(rbind, sim[26*(0:(N_simulations-1))+11])
  HIV_SW = do.call(rbind, sim[26*(0:(N_simulations-1))+12])
  HIV_MP = do.call(rbind, sim[26*(0:(N_simulations-1))+13])
  HIV_MSM = do.call(rbind, sim[26*(0:(N_simulations-1))+14])
  ART_GP = do.call(rbind, sim[26*(0:(N_simulations-1))+15])
  ART_SW = do.call(rbind, sim[26*(0:(N_simulations-1))+16])
  ART_MP = do.call(rbind, sim[26*(0:(N_simulations-1))+17])
  ART_MSM = do.call(rbind, sim[26*(0:(N_simulations-1))+18])
  HIV_GP_prop = do.call(rbind, sim[26*(0:(N_simulations-1))+19])
  HIV_SW_prop = do.call(rbind, sim[26*(0:(N_simulations-1))+20])
  HIV_MP_prop = do.call(rbind, sim[26*(0:(N_simulations-1))+21])
  HIV_MSM_prop = do.call(rbind, sim[26*(0:(N_simulations-1))+22])
  ART_GP_prop = do.call(rbind, sim[26*(0:(N_simulations-1))+23])
  ART_SW_prop = do.call(rbind, sim[26*(0:(N_simulations-1))+24])
  ART_MP_prop = do.call(rbind, sim[26*(0:(N_simulations-1))+25])
  ART_MSM_prop = do.call(rbind, sim[26*(0:(N_simulations-1))+26])
  
  N_years = floor(dim(HIV_counts)[2]/4)
  
  
  fetch_state_count = function(A, i){
    A[,i]
  }
  
  
  fetch_state_prop = function(A, i){
    A[,i]/apply(A, 1, sum)*100
  }
  
  
  analysis = function(u){
    quantile(u, c(alpha/2, .5, 1-alpha/2))
  }
  
  
  times = yyyy_mm_dd(start_year, N_years)
  
  state_count_analysis = function(i){
    GP = do.call(rbind, lapply(GPs, fetch_state_count, i=i))
    SW = do.call(rbind, lapply(SWs, fetch_state_count, i=i))
    MP = do.call(rbind, lapply(MPs, fetch_state_count, i=i))
    MSM = do.call(rbind, lapply(MSMs, fetch_state_count, i=i))
    
    GP_count_stats = as.data.frame(t(apply(GP, 2, analysis)))
    GP_count_stats$sub_pop = "GP"
    GP_count_stats$State = i
    GP_count_stats$Time = times
    
    SW_count_stats = as.data.frame(t(apply(SW, 2, analysis)))
    SW_count_stats$sub_pop = "FSW"
    SW_count_stats$State = i
    SW_count_stats$Time = times
    
    MP_count_stats = as.data.frame(t(apply(MP, 2, analysis)))
    MP_count_stats$sub_pop = "SC"
    MP_count_stats$State = i
    MP_count_stats$Time = times
    
    MSM_count_stats = as.data.frame(t(apply(MSM, 2, analysis)))
    MSM_count_stats$sub_pop = "MSM"
    MSM_count_stats$State = i
    MSM_count_stats$Time = times
    
    out = rbind(GP_count_stats, SW_count_stats, MP_count_stats, MSM_count_stats)
    names(out)[1:3] = c("LL", "Median", "UL")
    out
  }
  
  state_prop_analysis = function(i){
    GP = do.call(rbind, lapply(GPs, fetch_state_prop, i=i))
    SW = do.call(rbind, lapply(SWs, fetch_state_prop, i=i))
    MP = do.call(rbind, lapply(MPs, fetch_state_prop, i=i))
    MSM = do.call(rbind, lapply(MSMs, fetch_state_prop, i=i))
    
    GP_prop_stats = as.data.frame(t(apply(GP, 2, analysis)))
    GP_prop_stats$sub_pop = "GP"
    GP_prop_stats$State = i
    GP_prop_stats$Time = times
    
    SW_prop_stats = as.data.frame(t(apply(SW, 2, analysis)))
    SW_prop_stats$sub_pop = "FSW"
    SW_prop_stats$State = i
    SW_prop_stats$Time = times
    
    MP_prop_stats = as.data.frame(t(apply(MP, 2, analysis)))
    MP_prop_stats$sub_pop = "SC"
    MP_prop_stats$State = i
    MP_prop_stats$Time = times
    
    MSM_prop_stats = as.data.frame(t(apply(MSM, 2, analysis)))
    MSM_prop_stats$sub_pop = "MSM"
    MSM_prop_stats$State = i
    MSM_prop_stats$Time = times
    
    out = rbind(GP_prop_stats, SW_prop_stats, MP_prop_stats, MSM_prop_stats)
    names(out)[1:3] = c("LL", "Median", "UL")
    out
  }
  
  incidence_analysis = function(i){
    incidences = 1000*do.call(rbind, lapply(incidences, fetch_state_count, i=i))
    incidence_stats = as.data.frame(t(apply(incidences, 2, analysis)))
    if(i==1) incidence_stats$sub_pop = "GP"
    if(i==2) incidence_stats$sub_pop = "FSW"
    if(i==3) incidence_stats$sub_pop = "SC"
    if(i==4) incidence_stats$sub_pop = "MSM"
    incidence_stats$Time = times
    names(incidence_stats)[1:3] = c("LL", "Median", "UL")
    
    incidence_stats
  }
  
  count_stats = state_count_analysis(1)
  prop_stats = state_prop_analysis(1)
  incidence_stats = incidence_analysis(1)
  for(i in 2:5){
    count_stats = rbind(count_stats, state_count_analysis(i))
    prop_stats = rbind(prop_stats, state_prop_analysis(i))
    if(i<=4) incidence_stats = rbind(incidence_stats, incidence_analysis(i))
  }
  
  overall_incidences = as.data.frame(t(apply(overall_incidences, 2, analysis)))
  overall_incidences$Time = times
  names(overall_incidences)[1:3] = c("LL", "Median", "UL")
  
  HIV_counts = as.data.frame(t(apply(HIV_counts, 2, analysis)))
  HIV_counts$Time = times
  names(HIV_counts)[1:3] = c("LL", "Median", "UL")
  
  HIV_props = as.data.frame(t(apply(HIV_props, 2, analysis)))
  HIV_props$Time = times
  names(HIV_props)[1:3] = c("LL", "Median", "UL")
  
  ART_counts = as.data.frame(t(apply(ART_counts, 2, analysis)))
  ART_counts$Time = times
  names(ART_counts)[1:3] = c("LL", "Median", "UL")
  
  ART_props = as.data.frame(t(apply(ART_props, 2, analysis)))
  ART_props$Time = times
  names(ART_props)[1:3] = c("LL", "Median", "UL")
  
  HIV_counts_GP = as.data.frame(t(apply(HIV_GP, 2, analysis)))
  HIV_counts_GP$Time = times
  names(HIV_counts_GP)[1:3] = c("LL", "Median", "UL")
  HIV_counts_GP$sub_pop = "GP"
  
  HIV_counts_SW = as.data.frame(t(apply(HIV_SW, 2, analysis)))
  HIV_counts_SW$Time = times
  names(HIV_counts_SW)[1:3] = c("LL", "Median", "UL")
  HIV_counts_SW$sub_pop = "FSW"
  
  HIV_counts_MP = as.data.frame(t(apply(HIV_MP, 2, analysis)))
  HIV_counts_MP$Time = times
  names(HIV_counts_MP)[1:3] = c("LL", "Median", "UL")
  HIV_counts_MP$sub_pop = "SC"
  
  HIV_counts_MSM = as.data.frame(t(apply(HIV_MSM, 2, analysis)))
  HIV_counts_MSM$Time = times
  names(HIV_counts_MSM)[1:3] = c("LL", "Median", "UL")
  HIV_counts_MSM$sub_pop = "MSM"
  
  HIV_counts_subgrp = rbind(HIV_counts_GP, HIV_counts_SW, HIV_counts_MP, HIV_counts_MSM)
  
  HIV_props_GP = as.data.frame(t(apply(HIV_GP_prop, 2, analysis)))
  HIV_props_GP$Time = times
  names(HIV_props_GP)[1:3] = c("LL", "Median", "UL")
  HIV_props_GP$sub_pop = "GP"
  
  HIV_props_SW = as.data.frame(t(apply(HIV_SW_prop, 2, analysis)))
  HIV_props_SW$Time = times
  names(HIV_props_SW)[1:3] = c("LL", "Median", "UL")
  HIV_props_SW$sub_pop = "FSW"
  
  HIV_props_MP = as.data.frame(t(apply(HIV_MP_prop, 2, analysis)))
  HIV_props_MP$Time = times
  names(HIV_props_MP)[1:3] = c("LL", "Median", "UL")
  HIV_props_MP$sub_pop = "SC"
  
  HIV_props_MSM = as.data.frame(t(apply(HIV_MSM_prop, 2, analysis)))
  HIV_props_MSM$Time = times
  names(HIV_props_MSM)[1:3] = c("LL", "Median", "UL")
  HIV_props_MSM$sub_pop = "MSM"
  
  HIV_props_subgrp = rbind(HIV_props_GP, HIV_props_SW, HIV_props_MP, HIV_props_MSM)
  
  ART_counts_GP = as.data.frame(t(apply(ART_GP, 2, analysis)))
  ART_counts_GP$Time = times
  names(ART_counts_GP)[1:3] = c("LL", "Median", "UL")
  ART_counts_GP$sub_pop = "GP"
  
  ART_counts_SW = as.data.frame(t(apply(ART_SW, 2, analysis)))
  ART_counts_SW$Time = times
  names(ART_counts_SW)[1:3] = c("LL", "Median", "UL")
  ART_counts_SW$sub_pop = "FSW"
  
  ART_counts_MP = as.data.frame(t(apply(ART_MP, 2, analysis)))
  ART_counts_MP$Time = times
  names(ART_counts_MP)[1:3] = c("LL", "Median", "UL")
  ART_counts_MP$sub_pop = "SC"
  
  ART_counts_MSM = as.data.frame(t(apply(ART_MSM, 2, analysis)))
  ART_counts_MSM$Time = times
  names(ART_counts_MSM)[1:3] = c("LL", "Median", "UL")
  ART_counts_MSM$sub_pop = "MSM"
  
  ART_counts_subgrp = rbind(ART_counts_GP, ART_counts_SW, ART_counts_MP, ART_counts_MSM)
  
  ART_props_GP = as.data.frame(t(apply(ART_GP_prop, 2, analysis)))
  ART_props_GP$Time = times
  names(ART_props_GP)[1:3] = c("LL", "Median", "UL")
  ART_props_GP$sub_pop = "GP"
  
  ART_props_SW = as.data.frame(t(apply(ART_SW_prop, 2, analysis)))
  ART_props_SW$Time = times
  names(ART_props_SW)[1:3] = c("LL", "Median", "UL")
  ART_props_SW$sub_pop = "FSW"
  
  ART_props_MP = as.data.frame(t(apply(ART_MP_prop, 2, analysis)))
  ART_props_MP$Time = times
  names(ART_props_MP)[1:3] = c("LL", "Median", "UL")
  ART_props_MP$sub_pop = "SC"
  
  ART_props_MSM = as.data.frame(t(apply(ART_MSM_prop, 2, analysis)))
  ART_props_MSM$Time = times
  names(ART_props_MSM)[1:3] = c("LL", "Median", "UL")
  ART_props_MSM$sub_pop = "MSM"
  
  ART_props_subgrp = rbind(ART_props_GP, ART_props_SW, ART_props_MP, ART_props_MSM)
  
  return(list(count_stats=count_stats, prop_stats=prop_stats, 
              HIV_counts_all=HIV_counts, HIV_props_all=HIV_props,
              ART_counts_all=ART_counts, ART_props_all=ART_props,
              incidence_stats=incidence_stats, 
              overall_incidences=overall_incidences,
              HIV_counts=HIV_counts_subgrp,
              ART_counts=ART_counts_subgrp,
              HIV_props=HIV_props_subgrp,
              ART_props=ART_props_subgrp))
}





plot_outcomes_one = function(sim_results, stats_type, state, sub_pop){
  data = switch(stats_type,
                "props" = {sim_results$prop_stats},
                "incidences" = {sim_results$incidence_stats},
                "overall_incidence" = {sim_results$overall_incidences},
                "counts" = {sim_results$count_stats},
                "HIV_counts"= {sim_results$HIV_counts},
                "ART_counts"= {sim_results$ART_counts},
                "HIV_props"= {sim_results$HIV_props},
                "ART_props"= {sim_results$ART_props},
                "HIV_counts_all"= {sim_results$HIV_counts_all},
                "HIV_props_all"= {sim_results$HIV_props_all},
                "ART_counts_all"= {sim_results$ART_counts_all},
                "ART_props_all"= {sim_results$ART_props_all}
  )
  
  if(stats_type %in% c("overall_incidence", 
                       "HIV_counts_all", "HIV_props_all",
                       "ART_counts_all", "ART_props_all")){
    data = data
  } else if(stats_type %in% c("incidences", "HIV_counts", "ART_counts",
                              "HIV_props", "ART_props")){
    data = data[data$sub_pop==sub_pop,]
  } else{
    data = data[data$State==state & data$sub_pop==sub_pop,]
  }
  
  
  sub_pop_str = switch(sub_pop[1],
                       "FSW" = {": female sex workers"},
                       "SC" = {": male partners"},
                       "GP" = {": general population"},
                       "MSM" = {": MSM"}
  )
  
  state_str = switch(state,
                     "1" = {"HIV negative "},
                     "2" = {"Undiagnosed HIV "},
                     "3" = {"Diagnosed HIV Pre-ART "},
                     "4" = {"Viral load suppressed HIV on ART "},
                     "5" = {"ART Failed/Discontinued "}
  )
  
  
  stats_type_str = switch(stats_type,
                          "props" = {"Proportion of "},
                          "incidences" = {"Incidence "},
                          "overall_incidence" = {"Incidence "},
                          "counts" = {"Number of "},
                          "HIV_counts" = {"Number of people living with HIV "},
                          "HIV_counts_all" = {"Number of people living with HIV "},
                          "HIV_props" = {"Proportion of people living with HIV "},
                          "HIV_props_all" = {"Proportion of people living with HIV "},
                          "ART_counts" = {"Number of people on antiretroviral therapy "},
                          "ART_counts_all" = {"Number of people on antiretroviral therapy "},
                          "ART_props_all" = {"Proportion of people with HIV who are on ART "},
                          "ART_props" = {"Proportion of people with HIV who are on ART "}
  )
  
  stats_type_str2 = ifelse(stats_type == "counts", "by year (no. of cases)", "by year")
  
  
  if(stats_type %in% c("overall_incidence", 
                       "HIV_counts_all", "HIV_props_all",
                       "ART_counts_all", "ART_props_all")){ 
    str = paste0(stats_type_str, stats_type_str2)
  } else if(stats_type %in% c("incidences", "HIV_counts", "ART_counts",
                              "HIV_props", "ART_props")){
    str = paste0(stats_type_str, stats_type_str2, sub_pop_str)
  } else{
    str = paste0(stats_type_str, state_str, stats_type_str2, sub_pop_str)
  }
  
  pol_coords_x = c(data$Time, sort(data$Time, decreasing=T))
  pol_coords_y = c(data$LL, data$UL[length(data$UL):1])
  
  if(stats_type %in% c("counts", "HIV_counts", "HIV_counts_all", 
                       "ART_counts", "ART_counts_all")){
    ylab_str = "No. of cases"
  } else if(stats_type %in% c("incidences", "overall_incidence")){
    ylab_str = "No. of cases/1000 person years"
  } else{
    ylab_str = "%"
  }
  
  
  p = ggplot()  + 
    # scale_x_continuous(expand = expansion(c(0,0.025))) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          axis.text=element_text(size=12), 
          axis.title=element_text(size=14, face="bold"),
          axis.line=element_line(),
          plot.title = element_text(size=18, 
                                    hjust = 0.5, 
                                    face="bold")) +
    geom_polygon(aes(pol_coords_x, pol_coords_y), 
                 fill = "orange", colour = "skyblue", 
                 alpha = 0.5) +
    geom_line(aes(x=Time, y=Median), data=data, 
              colour="dodgerblue2", linewidth=1.5) + 
    geom_point(aes(Time, Median), data=data, size=3, 
               colour="#CC0000", shape=21, fill="white") + 
    ggtitle(str) + labs(x="Year", y=ylab_str) 
  
  p
}





plot_outcomes_some = function(sim_results, stats_type, state, sub_pops){
  data = switch(stats_type,
                "props" = {sim_results$prop_stats},
                "incidences" = {sim_results$incidence_stats},
                "HIV_counts"= {sim_results$HIV_counts},
                "ART_counts"= {sim_results$ART_counts},
                "HIV_props"= {sim_results$HIV_props},
                "ART_props"= {sim_results$ART_props},
                "overall_incidence" = {sim_results$overall_incidences},
                "counts" = {sim_results$count_stats}
  )
  
  if(stats_type %in% c("incidences", "HIV_counts", "ART_counts",
                       "HIV_props", "ART_props")){
    data = data[data$sub_pop %in% sub_pops,]
  } else{
    data = data[data$State==state & data$sub_pop %in% sub_pops,]
  }
  
  
  
  state_str = switch(state,
                     "1" = {"HIV negative"},
                     "2" = {"Undiagnosed HIV"},
                     "3" = {"Diagnosed HIV Pre-ART"},
                     "4" = {"Viral load suppressed HIV on ART"},
                     "5" = {"ART Failed/Discontinued"}
  )
  
  
  stats_type_str = switch(stats_type,
                          "props" = {"Proportion of "},
                          "incidences" = {"Incidence"},
                          "HIV_counts" = {"Number of people living with HIV"},
                          "HIV_props" = {"Proportion of people living with HIV"},
                          "ART_counts" = {"Number of people on antiretroviral therapy"},
                          "ART_props" = {"Proportion of people with HIV who are on ART"},
                          "counts" = {""}
  )
  
  stats_type_str2 = ifelse(stats_type == "counts", " by year (no. of cases)", " by year")
  
  
  
  if(stats_type %in% c("incidences", "HIV_counts", "ART_counts",
                       "HIV_props", "ART_props")){
    str = paste0(stats_type_str,stats_type_str2)
  } else{
    str = paste0(stats_type_str, state_str, stats_type_str2)
  }
  
  
  if(stats_type %in% c("counts", "HIV_counts", "ART_counts")){
    ylab_str = "No. of cases"
  } else if(stats_type %in% c("incidences", "overall_incidence")){
    ylab_str = "No. of cases/1000 person years"
  } else{
    ylab_str = "%"
  }
  
  #data$sub_pop = gsub("MP", "SC", data$sub_pop)
  data$Time = as.Date(data$Time)
  
  p = ggplot(data, aes(x=Time, y=Median, group=sub_pop, color=sub_pop)) + 
    geom_line(linewidth=1.5) +
    geom_errorbar(aes(ymin=LL, ymax=UL), width=.2, linewidth=1.5,
                  position=position_dodge(0.05), colour="black") +
    geom_point(size=3, shape=21, fill="white") +
    labs(col="Group") + 
    ggtitle(str) + labs(x="Year", y=ylab_str) + 
    # scale_x_continuous(expand = expansion(c(0,0.025))) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          axis.text=element_text(size=10), 
          axis.title=element_text(size=12, face="bold"),
          axis.line=element_line(),
          plot.title = element_text(size=14, hjust = 0.5, face="bold")) 
  
  p
}





plot_outcomes = function(sim_results, stats_type, state, sub_pops){
  if(length(sub_pops)==1 | stats_type %in% c("overall_incidence", 
                                             "HIV_counts_all", "HIV_props_all",
                                             "ART_counts_all", "ART_props_all")){
    plot_outcomes_one(sim_results, stats_type, state, sub_pops)
  } else{
    plot_outcomes_some(sim_results, stats_type, state, sub_pops)
  }
}





group_lumping = function(sub_pops, lump, sim_results, stats_type, state){
  metric_data = switch(stats_type,
                       "props" = {sim_results$prop_stats},
                       "incidences" = {sim_results$incidence_stats},
                       "HIV_counts"= {sim_results$HIV_counts},
                       "ART_counts"= {sim_results$ART_counts},
                       "HIV_props"= {sim_results$HIV_props},
                       "ART_props"= {sim_results$ART_props},
                       "overall_incidence" = {sim_results$overall_incidences},
                       "counts" = {sim_results$count_stats},
                       "HIV_counts_all"= {sim_results$HIV_counts_all},
                       "HIV_props_all"= {sim_results$HIV_props_all},
                       "ART_counts_all"= {sim_results$ART_counts_all},
                       "ART_props_all"= {sim_results$ART_props_all}
  )
  
  if(!(stats_type %in% c("HIV_counts_all", "HIV_props_all",
                         "ART_counts_all", "ART_props_all",
                         "overall_incidence"))){
    metric_data = metric_data[metric_data$sub_pop %in% sub_pops, ]
    lump = lump[lump %in% sub_pops]
  }
  
  if(stats_type %in% c("incidences", "HIV_counts", "ART_counts",
                       "HIV_props", "ART_props")){
    leftovers = metric_data[!(metric_data$sub_pop %in% lump), -c(1,3)]
    metric_data = metric_data[metric_data$sub_pop %in% lump, -c(1,3)]
  } else if(!(stats_type %in% c("HIV_counts_all", "HIV_props_all",
                                "ART_counts_all", "ART_props_all",
                                "overall_incidence"))){
    leftovers = metric_data[metric_data$State==state & 
                              !(metric_data$sub_pop %in% lump), -c(1,3)]
    metric_data = metric_data[metric_data$State==state &
                                metric_data$sub_pop %in% lump, -c(1,3)]
  } else{
    metric_data = metric_data[, c(4,2)]
  }
  
  if(!(stats_type %in% c("HIV_counts_all", "HIV_props_all",
                         "ART_counts_all", "ART_props_all",
                         "overall_incidence"))){
    counts_data = sim_results$count_stats[sim_results$count_stats$sub_pop %in% lump, -c(1,3)] %>%
      dplyr::group_by(sub_pop, Time) %>%
      summarize(N = sum(Median, na.rm = TRUE)) %>%
      as.data.frame
    
    Totals = counts_data %>% 
      dplyr::group_by(Time) %>%
      summarize(N = sum(N, na.rm = TRUE)) %>%
      as.data.frame
    
    Totals1 = Totals
    if(length(lump) >= 2){
      for(i in 1:(length(lump) - 1)) Totals1 = rbind(Totals1, Totals)
    }
    
    weights = counts_data$N/Totals1$N
    metric_data$Weighted = metric_data$Median*weights
    
    df = metric_data %>% 
      dplyr::group_by(Time) %>%
      summarize(Value = sum(Weighted, na.rm = TRUE)) %>%
      as.data.frame
    
    if(length(lump) >= 2){
      if(length(lump) == length(sub_pops) && prod(sort(lump) == sort(sub_pops))){
        df$sub_pop = paste(lump, collapse = " + ")
      } else{
        df$sub_pop = "Other"
      }
    } else if(length(lump) == 0){
      df$sub_pop = character(0)
    } else{
      df$sub_pop = lump
    }
    
    df2 = leftovers[,c("Time", "Median", "sub_pop")]
    names(df2) = names(df)
    df = rbind(df, df2)
  } else{
    df = metric_data
  }
  
  df
}







lumped_groups_plot = function(sub_pops, lump, sim_results, stats_type, state){
  df = group_lumping(sub_pops, lump, sim_results, stats_type, state)
  
  state_str = switch(state,
                     "1" = {"HIV negative"},
                     "2" = {"Undiagnosed HIV"},
                     "3" = {"Diagnosed HIV Pre-ART"},
                     "4" = {"Viral load suppressed HIV on ART"},
                     "5" = {"ART Failed/Discontinued"}
  )
  
  
  stats_type_str = switch(stats_type,
                          "props" = {"Proportion of "},
                          "incidences" = {"Incidence"},
                          "HIV_counts" = {"Number of people living with HIV"},
                          "HIV_props" = {"Proportion of people living with HIV"},
                          "ART_counts" = {"Number of people on antiretroviral therapy"},
                          "ART_props" = {"Proportion of people with HIV who are on ART"},
                          "counts" = {""}
  )
  
  stats_type_str2 = ifelse(stats_type == "counts", " by year (no. of cases)", " by year")
  
  if(stats_type %in% c("incidences", "HIV_counts", "ART_counts",
                       "HIV_props", "ART_props")){
    str = paste0(stats_type_str,stats_type_str2)
  } else{
    str = paste0(stats_type_str, state_str, stats_type_str2)
  }
  
  
  if(stats_type %in% c("counts", "HIV_counts", "ART_counts")){
    ylab_str = "No. of cases"
  } else if(stats_type %in% c("incidences", "overall_incidence")){
    ylab_str = "No. of cases/1000 person years"
  } else{
    ylab_str = "%"
  }
  
  #df$sub_pop = gsub("MP", "SC", df$sub_pop)
  subpop = as.factor(df$sub_pop)
  lev = levels(subpop)
  if("Other" %in% lev){
    temp = lev[lev != "Other"]
    subpop = factor(df$sub_pop, levels = c(temp, "Other"))
  }
  df$sub_pop = subpop
  
  df$Time = as.Date(df$Time)
  
  p = ggplot(df, aes(x=Time, y=Value, group=sub_pop, color=sub_pop)) + 
    geom_line(linewidth=1.5) +
    geom_point(size=3, shape=21, fill="white") +
    labs(col="Group") +
    ggtitle(str) + 
    labs(x="Year", y=ylab_str) + 
    # scale_x_continuous(expand = expansion(c(0,0.025))) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          axis.text=element_text(size=12), 
          axis.title=element_text(size=14, face="bold"),
          axis.line=element_line(),
          legend.text=element_text(size=12),
          legend.title=element_text(size=13, face="bold"),
          plot.title = element_text(size=18, 
                                    hjust = 0.5, 
                                    face="bold")) #+ 
  #expand_limits(y = 0)
  
  p
}






With_Without_PrEP_plot = function(sub_pops, sim_results_PrEP, 
                                  sim_results_No_PrEP,  
                                  stats_type, state){
  stats_fetch = function(sim_results, stats_type){
    switch(stats_type,
           "props" = {sim_results$prop_stats},
           "incidences" = {sim_results$incidence_stats},
           "HIV_counts"= {sim_results$HIV_counts},
           "ART_counts"= {sim_results$ART_counts},
           "HIV_props"= {sim_results$HIV_props},
           "ART_props"= {sim_results$ART_props},
           "overall_incidence" = {sim_results$overall_incidences},
           "counts" = {sim_results$count_stats},
           "HIV_counts_all"= {sim_results$HIV_counts_all},
           "HIV_props_all"= {sim_results$HIV_props_all},
           "ART_counts_all"= {sim_results$ART_counts_all},
           "ART_props_all"= {sim_results$ART_props_all}
    )
  }
  
  stats_PrEP = stats_fetch(sim_results_PrEP, stats_type)
  stats_No_PrEP = stats_fetch(sim_results_No_PrEP, stats_type)
  if(!(stats_type %in% c("overall_incidence", "HIV_counts_all", 
                         "HIV_props_all", "ART_counts_all", "ART_props_all"))){
    stats_PrEP = stats_PrEP[stats_PrEP$sub_pop %in% sub_pops,]
    stats_No_PrEP = stats_No_PrEP[stats_No_PrEP$sub_pop %in% sub_pops,]
  }
  
  if(stats_type %in% c("props", "counts")){
    stats_PrEP = stats_PrEP[stats_PrEP$State == state, ]
    stats_No_PrEP = stats_No_PrEP[stats_No_PrEP$State == state, ]
  }
  
  stats_PrEP$sub_pop = paste(stats_PrEP$sub_pop, "PrEP")
  stats_No_PrEP$sub_pop = paste(stats_No_PrEP$sub_pop, "No PrEP")
  stats_all = rbind(stats_No_PrEP, stats_PrEP)
  stats_all$cols = unlist(lapply(strsplit(stats_all$sub_pop, " "), function(x) x[1]))
  stats_all$ltypes = unlist(lapply(strsplit(stats_all$sub_pop, " "), function(x) x[2]))
  
  stats_all$ltypes[stats_all$ltypes == "PrEP"] = 0
  stats_all$ltypes[stats_all$ltypes != 0] = 1
  data = stats_all 
  # data[,c("sub_pop", "cols")] = apply(data[,c("sub_pop", "cols")], 
  #                                     2, function(x){gsub("MP", "SC", x)})
  data$Time = as.Date(data$Time)
  data$Median = as.numeric(data$Median)
  
  state_str = switch(state,
                     "1" = {"HIV negative"},
                     "2" = {"Undiagnosed HIV"},
                     "3" = {"Diagnosed HIV Pre-ART"},
                     "4" = {"Viral load suppressed HIV on ART"},
                     "5" = {"ART Failed/Discontinued"}
  )
  
  
  stats_type_str = switch(stats_type,
                          "props" = {"Proportion of "},
                          "incidences" = {"Incidence"},
                          "HIV_counts" = {"Number of people living with HIV"},
                          "HIV_props" = {"Proportion of people living with HIV"},
                          "ART_counts" = {"Number of people on antiretroviral therapy"},
                          "ART_props" = {"Proportion of people with HIV who are on ART"},
                          "counts" = {""},
                          "HIV_props_all" = {"Proportion of people living with HIV"},
                          "overall_incidence" = {"Incidence"},
                          "HIV_counts_all" = {"Number of people living with HIV "},
                          "ART_counts_all" = {"Number of people on antiretroviral therapy "},
                          "ART_props_all" = {"Proportion of people with HIV who are on ART "}
  )
  
  stats_type_str2 = ifelse(stats_type == "counts", " by year (no. of cases)", " by year")
  
  if(stats_type %in% c("incidences", "HIV_counts", "ART_counts",
                       "HIV_props", "ART_props", "overall_HIV_count",
                       "HIV_props_all", "overall_incidence", 
                       "HIV_counts_all", "ART_counts_all", "ART_props_all")){
    str = paste0(stats_type_str, stats_type_str2)
  } else{
    str = paste0(stats_type_str, state_str, stats_type_str2)
  }
  
  
  if(stats_type %in% c("counts", "HIV_counts", "ART_counts", "HIV_counts_all")){
    ylab_str = "No. of cases"
  } else if(stats_type %in% c("incidences", "overall_incidence")){
    ylab_str = "No. of cases/1000 person years"
  } else{
    ylab_str = "%"
  }
  
  
  if(!(stats_type %in% c("HIV_props_all", 
                         "overall_incidence",
                         "HIV_counts_all",
                         "ART_props_all",
                         "ART_counts_all"))){
    p = ggplot(data, aes(x=Time, y=Median, group=sub_pop, color=cols))
  } else{
    p = ggplot(data, aes(x=Time, y=Median))
  }
  p = p + 
    geom_line(linewidth=1.5, 
              aes(linetype = factor(ltypes, labels=c("Yes", "No")))) +
    geom_point(size=1, shape=21, fill="white") +
    labs(col="Group", linetype ="PrEP") + 
    ggtitle(str) + 
    # scale_x_continuous(expand = expansion(c(0,0.025))) +
    labs(x="Year", y=ylab_str) +
    theme(panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          axis.text=element_text(size=12), 
          axis.title=element_text(size=14, face="bold"),
          axis.line=element_line(),
          legend.text=element_text(size=12),
          legend.title=element_text(size=13, face="bold"),
          plot.title = element_text(size=18, 
                                    hjust = 0.5, 
                                    face="bold"))
  
  
  p
}






With_Without_PrEP_lumped_plot = function(sub_pops, lump, sim_results_PrEP, 
                                         sim_results_No_PrEP, stats_type, state){
  stats_PrEP = group_lumping(sub_pops, lump, sim_results_PrEP, stats_type, state)
  stats_No_PrEP = group_lumping(sub_pops, lump, sim_results_No_PrEP, stats_type, state)
  stats_PrEP$sub_pop = paste(stats_PrEP$sub_pop, "PrEP")
  stats_No_PrEP$sub_pop = paste(stats_No_PrEP$sub_pop, "No PrEP")
  stats_all = rbind(stats_No_PrEP, stats_PrEP)
  
  stats_all$cols = unlist(lapply(strsplit(stats_all$sub_pop, " "), function(x) x[1]))
  stats_all$ltypes = unlist(lapply(strsplit(stats_all$sub_pop, " "), function(x) x[2]))
  stats_all$ltypes[stats_all$ltypes == "PrEP"] = 0
  stats_all$ltypes[stats_all$ltypes != 0] = 1
  data = stats_all 
  #data[,3] = gsub("MP", "SC", data[,3])
  data$Time = as.Date(data$Time)
  
  state_str = switch(state,
                     "1" = {"HIV negative"},
                     "2" = {"Undiagnosed HIV"},
                     "3" = {"Diagnosed HIV Pre-ART"},
                     "4" = {"Viral load suppressed HIV on ART"},
                     "5" = {"ART Failed/Discontinued"}
  )
  
  
  stats_type_str = switch(stats_type,
                          "props" = {"Proportion of "},
                          "incidences" = {"Incidence"},
                          "overall_incidences" = {"Incidence"},
                          "HIV_counts" = {"Number of people living with HIV"},
                          "HIV_counts_all" = {"Number of people living with HIV"},
                          "HIV_props" = {"Proportion of people living with HIV"},
                          "ART_counts" = {"Number of people on antiretroviral therapy"},
                          "ART_props" = {"Proportion of people with HIV who are on ART"},
                          "counts" = {""},
                          "HIV_props_all" = {"Proportion of people living with HIV "},
                          "ART_counts_all" = {"Number of people on antiretroviral therapy "},
                          "ART_props_all" = {"Proportion of people with HIV who are on ART "}
  )
  
  stats_type_str2 = ifelse(stats_type == "counts", " by year (no. of cases)", " by year")
  
  
  
  if(stats_type %in% c("incidences", "HIV_counts", "ART_counts",
                       "HIV_props", "ART_props")){
    str = paste0(stats_type_str,stats_type_str2)
  } else{
    str = paste0(stats_type_str, state_str, stats_type_str2)
  }
  
  
  if(stats_type %in% c("counts", "HIV_counts", "ART_counts")){
    ylab_str = "No. of cases"
  } else if(stats_type %in% c("incidences", "overall_incidence")){
    ylab_str = "No. of cases/1000 person years"
  } else{
    ylab_str = "%"
  }
  
  p = ggplot(data, aes(x=Time, y=Value, group=sub_pop, color=cols)) + 
    geom_line(linewidth=1.5, 
              aes(linetype = factor(ltypes, labels=c("Yes", "No")))) +
    geom_point(size=1, shape=21, fill="white") +
    # scale_x_continuous(expand = expansion(c(0,0.025))) + 
    labs(col="Group", linetype ="PrEP") + 
    ggtitle(str) + 
    labs(x="Year", y=ylab_str) + 
    theme(panel.background = element_blank(),
          panel.grid.major = element_line(colour = "grey92"),
          axis.text=element_text(size=12), 
          axis.title=element_text(size=14, face="bold"),
          legend.text=element_text(size=12),
          legend.title=element_text(size=13, face="bold"),
          axis.line=element_line(),
          plot.title = element_text(size=18, 
                                    hjust = 0.5, 
                                    face="bold")) 
  
  
  p
}







populate = function(pop_size, p_SW, p_MP, p_hiv_overall,
                    p_hiv_SW, p_hiv_MP, p_test_SW,
                    p_test_MP, p_test_overall,
                    p_ART_SW, p_ART_MP, p_ART_overall,
                    p_ART_sup, p_ART_adhere, p_MSM,
                    p_hiv_MSM, p_test_MSM, p_ART_MSM){
  SW = pop_size*p_SW
  MP = pop_size*p_MP
  MSM = pop_size*p_MSM
  GP = pop_size*(1 - p_SW - p_MP - p_MSM)
  
  hiv_SW = SW*p_hiv_SW
  hiv_MP = MP*p_hiv_MP
  hiv_MSM = MSM*p_hiv_MSM
  hiv_overall = pop_size*p_hiv_overall
  hiv_GP = hiv_overall - hiv_MP - hiv_SW - hiv_MSM
  
  hiv_neg_SW = SW - hiv_SW
  hiv_neg_MP = MP - hiv_MP
  hiv_neg_GP = GP - hiv_GP
  hiv_neg_MSM = MSM - hiv_MSM
  
  hiv_undiag_SW = hiv_SW*(1-p_test_SW)
  hiv_undiag_MP = hiv_MP*(1-p_test_MP)
  hiv_undiag_MSM = hiv_MSM*(1-p_test_MSM)
  hiv_undiag_overall = hiv_overall*(1-p_test_overall)
  hiv_undiag_GP = hiv_undiag_overall - hiv_undiag_MP - hiv_undiag_SW - hiv_undiag_MSM
  
  hiv_diag_SW = hiv_SW - hiv_undiag_SW
  hiv_diag_MP = hiv_MP - hiv_undiag_MP
  hiv_diag_MSM = hiv_MSM - hiv_undiag_MSM
  hiv_diag_overall = hiv_overall - hiv_undiag_overall
  
  hiv_ART_SW = hiv_diag_SW*p_ART_SW
  hiv_ART_MP = hiv_diag_MP*p_ART_MP
  hiv_ART_MSM = hiv_diag_MSM*p_ART_MSM
  hiv_ART_overall = hiv_diag_overall*p_ART_overall
  
  ART_sup_SW = hiv_ART_SW*p_ART_sup*p_ART_adhere
  ART_sup_MP = hiv_ART_MP*p_ART_sup*p_ART_adhere
  ART_sup_MSM = hiv_ART_MSM*p_ART_sup*p_ART_adhere
  ART_sup_overall = hiv_ART_overall*p_ART_sup*p_ART_adhere
  ART_sup_GP = ART_sup_overall - ART_sup_SW - ART_sup_MP - ART_sup_MSM
  
  p_ART_no_sup = 1 - p_ART_sup*p_ART_adhere
  ART_no_sup_SW = hiv_ART_SW*p_ART_no_sup
  ART_no_sup_MP = hiv_ART_MP*p_ART_no_sup
  ART_no_sup_MSM = hiv_ART_MSM*p_ART_no_sup
  ART_no_sup_overall = hiv_ART_overall*p_ART_no_sup
  ART_no_sup_GP = ART_no_sup_overall - ART_no_sup_SW - ART_no_sup_MP - ART_no_sup_MSM
  
  No_ART_SW = SW - hiv_neg_SW - hiv_undiag_SW - ART_sup_SW - ART_no_sup_SW
  No_ART_MP = MP - hiv_neg_MP - hiv_undiag_MP - ART_sup_MP - ART_no_sup_MP
  No_ART_MSM = MSM - hiv_neg_MSM - hiv_undiag_MSM - ART_sup_MSM - ART_no_sup_MSM
  No_ART_GP = GP - hiv_neg_GP - hiv_undiag_GP - ART_sup_GP - ART_no_sup_GP - No_ART_MSM
  
  SW_counts = round(c(hiv_neg_SW, hiv_undiag_SW, No_ART_SW, ART_sup_SW, ART_no_sup_SW))
  MP_counts = round(c(hiv_neg_MP, hiv_undiag_MP, No_ART_MP, ART_sup_MP, ART_no_sup_MP))
  MSM_counts = round(c(hiv_neg_MSM, hiv_undiag_MSM, No_ART_MSM, ART_sup_MSM, ART_no_sup_MSM))
  GP_counts = round(c(hiv_neg_GP, hiv_undiag_GP, No_ART_GP, ART_sup_GP, ART_no_sup_GP))
  
  return(list(SW_counts=SW_counts, MP_counts=MP_counts, 
              GP_counts=GP_counts, MSM_counts=MSM_counts))
}






# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots = c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout = matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}