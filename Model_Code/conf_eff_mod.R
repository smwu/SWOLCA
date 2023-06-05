convert_from_reference <- function(beta) {
  # Convert from b1 + b2*I(C=2) + b3*I(C=3) + b4*I(S=2) + b5*I(C=2,S=2) + b6*I(C=3,S=2)
  # to xi1*I(S=1,C=1) + xi2*I(S=1,C=2) + xi3*I(S=1,C=3) + xi4*I(S=2,C=1) + xi5*I(S=2,C=2) + xi6*I(S=2,C=3)
  xi_vec <- numeric(length(beta))
  xi_vec[1] <- xi_vec[1]
  xi_vec[2] <- beta[1] + beta[2]
  xi_vec[3] <- beta[1] + beta[3]
  xi_vec[4] <- beta[1] + beta[4]
  xi_vec[5] <- beta[1] + beta[2] + beta[4] + beta[5]
  xi_vec[6] <- beta[1] + beta[3] + beta[4] + beta[6]
  xi <- matrix(xi_vec, nrow = 3, ncol = 2)
  return(xi)
}

convert_to_reference <- function(xi) {
  # Convert to b1 + b2*I(C=2) + b3*I(C=3) + b4*I(S=2) + b5*I(C=2,S=2) + b6*I(C=3,S=2)
  xi_vec <- c(xi)
  beta <- numeric(length(xi_vec))
  beta[1] <- xi_vec[1]
  beta[2] <- xi_vec[2] - xi_vec[1]
  beta[3] <- xi_vec[3] - xi_vec[1]
  beta[4] <- xi_vec[4] - xi_vec[1]
  beta[5] <- xi_vec[5] - xi_vec[4] - beta[2]
  beta[6] <- xi_vec[6] - xi_vec[4] - beta[3]
  return(beta)
}

convert_to_reference_marg <- function(xi) {
  # Convert to b1 + b2*I(C=2) + b3*I(C=3)
  xi_vec <- c(xi)
  beta <- numeric(length(xi_vec))
  beta[1] <- xi_vec[1]
  beta[2] <- xi_vec[2] - xi_vec[1]
  beta[3] <- xi_vec[3] - xi_vec[1]
  return(beta)
}

glm_data <- data.frame(sim_pop$true_Ci, sim_pop$true_Si, sim_pop$Y_data)
colnames(glm_data) <- c("C", "S", "Y")
glm_data$C <- as.factor(glm_data$C)
glm_data$S <- as.factor(glm_data$S)
convert_to_reference(sim_pop$true_xi)
glm1 <- glm(Y ~ C * S, data = glm_data, family = binomial(link = "probit"))
summary(glm1)


marg_Phi <- get_marg_eff(sim_pop$true_xi, sim_pop$true_pi_s, sim_pop$true_pi, 
             sim_pop$N_s, sim_pop$N)
marg_xi <- qnorm(marg_Phi)
convert_to_reference_marg(marg_xi)
glm2 <- glm(Y ~ C, data = glm_data, family = binomial(link = "probit"))
summary(glm2)

glm3 <- glm(Y ~ C + S, data = glm_data, family = binomial(link = "probit"))
summary(glm3)

temp <- table(glm_data$C, glm_data$S)
marg_S <- colSums(temp)
marg_C <- rowSums(temp)
outer(marg_C, marg_S)/sum(marg_S)
temp
chisq.test(temp)

chisq.test(table(glm_data$C, glm_data$S))
chisq.test(table(glm_data$C, glm_data$Y))
chisq.test(table(glm_data$S, glm_data$Y))



# Simulate confounder
library(fastDummies)
s2 <- ifelse(glm_data$S == 2, 1, 0)
c2 <- ifelse(glm_data$C == 2, 1, 0)
c3 <- ifelse(glm_data$C == 3, 1, 0)
design_mat <- as.matrix(cbind(rep(1, 80000),c2, c3, s2))
beta_conf <- c(1.3, -1.1, -1.79, -0.82)
Phi_all <- pnorm(design_mat %*% as.matrix(beta_conf))
xi_conf <- convert_from_reference(c(beta_conf, 0, 0))
Y_conf <- rbinom(80000, 1, Phi_all)
glm_data$Y_conf <- Y_conf
summary(glm(Y_conf ~ C + S, data = glm_data, 
            family = binomial(link = "probit")))
chisq.test(table(glm_data$S, glm_data$Y_conf))

summary(glm(Y_conf ~ C, data = glm_data, 
            family = binomial(link = "probit")))
marg_Phi_conf <- get_marg_eff(xi_conf, sim_pop$true_pi_s, sim_pop$true_pi, 
                              sim_pop$N_s, sim_pop$N)
marg_Phi_conf
marg_Phi_conf_obs <- c(mean(glm_data$Y_conf[glm_data$C == 1]), 
                       mean(glm_data$Y_conf[glm_data$C == 2]), 
                       mean(glm_data$Y_conf[glm_data$C == 3]))
marg_xi_conf <- qnorm(marg_Phi)  # xi1*I(C=1) + xi2*I(C=2) + xi3*I(C=3)
convert_to_reference_marg(marg_xi_conf)
# compare
beta_conf

xi_under <- matrix(c(1, 0.3, -0.5, 0.5, -0.7, -1.3), nrow = 3, byrow = FALSE)
# Matrix of true outcome probabilities
true_Phi_mat <- pnorm(xi_under) 
# True outcome probabilities for each individual
true_Phi <- numeric(N)   
# Outcome data
Y_data <- numeric(N) 
for (s in 1:S) {
  for (k in 1:K) {
    # Individuals in stratum s class k
    N_s_k <- pop_inds[true_Si == s & true_Ci == k]
    with_outcome <- sample(N_s_k, round(length(N_s_k) * true_Phi_mat[k, s]))
    Y_data[with_outcome] <- 1
    true_Phi[N_s_k] <- true_Phi_mat[k, s]
  }
}
