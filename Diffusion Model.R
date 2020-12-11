# Ratcliff & Rouder's (1998) Diffusion Model

library("rtdists")
library(data.table)
library(lattice)
library(purrr)
library(latticeExtra)
library(dplyr)
library(tidyr)

data(rr98)
head(rr98)
data = rr98[!rr98$outlier == TRUE,] #removes outliers

# Takes the parameter values from Table 1 for Experiment 1 and creates a random distribution of response times

x<- seq(0, 1, length = 33)
probability_33 = (dnorm(x, 0.625, sd = 0.1875))/33
probability = (dnorm(x, 0.625, sd = 0.1875))


##
data_2 <- data  %>% group_by(id, instruction, strength) %>% 
  summarise(prop = mean(response == "light"), mean_rt = mean(rt), median_rt = mean(rt), sd_rt = 2*(sd(rt))) %>% 
  ungroup()

xyplot(prop ~ strength|id, data_2, group = instruction, type = "b", 
       auto.key = list(lines = TRUE), xlab = "Stimulus Value (Strength)",
       ylab = "Proportion of 'light' Responses"
         )
## Fig. 4
#z = mean(data$response == "light" & data$response_num == 1)
data_1 <- data  %>% group_by(id, instruction, strength) %>% 
  summarise(prop = mean(response_num == 1), mean_rt = 1000*(mean(rt)), 
            median_rt = mean(rt), sd_rt = 2*(sd(rt))) %>% 
  ungroup()

xyplot(mean_rt ~prop|id, data_1, group = instruction, type = "b", 
       auto.key = list(lines = FALSE), col = c("darkorchid1", "deepskyblue2"),
       xlab = "Response Probability",ylab = "Mean Reaction Time (ms)")

###### Fig 3B - Probability of each pixel proportion

plot(x, (dnorm(x, 0.375, sd = 0.1875))/33, type = "l", col = "darkorchid1", xlab = "Proportion of White Pixels", ylab = "Probability")
for (i in 1:1){lines(x,((dnorm(x, 0.625, sd = 0.1875))/33), col = "deepskyblue2")}
plot(x, pnorm(x, mean = 0.375, sd = 0.1875, lower.tail = TRUE, log.p = FALSE))

# do all people have the same drift rate? Isn't the probability distriution the same across all participants?
v = ((0.9*probability) - 0.45)/33 #divides probability by 33 
#v = (0.9*(dnorm(x, 0.625, sd = 0.1875))) - 0.45
plot(v)

xyplot(v ~ strength|id, data_2, group = instruction, type = "b", 
       auto.key = list(lines = TRUE), xlab = "Drift Rate",ylab = "Proportion of White Pixels")
       

########## Quantiles
quantiles <- c(0.1, 0.3, 0.5, 0.7, 0.9)
quantiles_data <- data %>%
  group_by(id, instruction, strength) %>%
  nest() %>% 
  mutate(quantiles = map(data, ~ as.data.frame(t(quantile(.x$rt, probs = quantiles))))) %>% 
  unnest(quantiles) %>% 
  gather("quantile", "rt",`10%`:`90%`) %>% 
  arrange(id, instruction, strength)

xyplot(rt*1000 ~ strength|id + instruction, quantiles_data, group = quantile, type = "b", 
       auto.key = list(lines = FALSE, space = "right"), ylab = "Mean Reaction Time (ms)", subset = instruction == "speed")

xyplot(rt*1000 ~ strength|id + instruction, quantiles_data, group = quantile, type = "b", 
       auto.key = list(lines = FALSE, space = "right"), ylab = "Mean Reaction Time (ms)", subset = instruction == "accuracy")

q_diff = qdiffusion(quantiles, a=0.079, v, z = 0.0395, t0=0.26, sv = 0.063, s = 0.1, response="u", scale_p = TRUE)
plot(quantiles, q_diff)
######## Modeling each participant?
# dataset to be fitted
data_fitted <- data %>% 
  group_by(id, instruction) %>% # we loop across both, id and instruction
  nest()
data_fitted
# function to assign drift rates to strengths 
objective_diffusion_separate <- function(pars, rt, response, drift, ...) {
  non_v_pars <- grep("^v", names(pars), invert = TRUE, value = TRUE)
  base_par <- length(non_v_pars)  # number of non-drift parameters
  densities <- vector("numeric", length(rt))
  for (i in seq_along(levels(drift))) {
    densities[drift == levels(drift)[i]] <- 
      ddiffusion(rt[drift == levels(drift)[i]], response=response[drift == levels(drift)[i]], 
                 a=pars["a"], t0=pars["t0"],  
                 sv=pars["sv"],
                 sz=if ("sz" %in% non_v_pars) pars["sz"] else 0.1,
                 z=if ("z" %in% non_v_pars) pars["z"]*pars["a"] else 0.5*pars["a"],
                 st0=if ("st0" %in% non_v_pars) pars["st0"] else 0, 
                 v=pars[base_par+i])
  }
  if (any(densities == 0)) return(1e6)
  return(-sum(log(densities)))
}

# function that creates random start values 
get_start <- function(base_par, n_drift = 5) {
  start1 <- c(
    a = runif(1, 0.5, 3),
    a_1 = runif(1, 0.5, 3), 
    a_2 = runif(1, 0.5, 3),
    t0 = runif(1, 0, 0.5), 
    z = runif(1, 0.4, 0.6), 
    sz = runif(1, 0, 0.5),
    sv = runif(1, 0, 0.5),
    st0 = runif(1, 0, 0.5),
    d = rnorm(1, 0, 0.05)
  )
  start2 <- sort(rnorm(n_drift), decreasing = FALSE)
  names(start2) <- paste0("v_", seq_len(n_drift))
  c(start1[base_par], start2)
}

# function that tries different random starting values until it works
ensure_fit <- 
  function(data, start_function, objective_function, 
           base_pars, n_drift = 5, n_fits = 1, 
           lower = c(rep(0, length(base_pars)), -Inf,
                     rep(-Inf,length(start_function(base_pars))-length(base_pars)))) {
    best_fit <- list(objective = 1e+06)
    for (i in seq_len(n_fits)) {
      start_ll <- 1e+06
      #browser()
      while(start_ll == 1e+06) {
        start <- start_function(base_pars, n_drift=n_drift)
        start_ll <- objective_function(start, 
                                       rt = data$rt, response = data$response_num, 
                                       drift = factor(data$strength, seq_len(n_drift)), 
                                       instruction = data$instruction)
      }
      cat("\nstart fitting.\n") 
      
      fit <- nlminb(start, objective_function, 
                    rt = data$rt, response = data$response_num, 
                    drift = factor(data, seq_len(n_drift)), 
                    instruction = data$instruction,
                    lower == "lower")
      
      if (fit$objective < best_fit$objective) best_fit <- fit
    }
    out <- as.data.frame(t(unlist(best_fit[1:3])))
    colnames(out) <- sub("par.", "", colnames(out))
    out
  }

fit_diffusion <- data_fitted %>% 
  mutate(fit = 
           map(data, 
               ~ensure_fit(data = ., start_function = get_start, 
                           objective_function = objective_diffusion_separate, 
                           as.numeric(response) %in% 1:2!!,
                            base_pars = c("a", "t0", "sv", "sz", "z")))) %>% 
  unnest(fit)

######
n_NH_acc = data[data$id == "nh" & data$instruction == "accuracy",] # number of accuracy trials
plot(n_NH_acc$rt~n_NH_acc$strength)
n_NH_speed = data[data$id == "nh" & data$instruction == "speed",] # number of speed trials

NH_speed_model = rdiffusion(4345, a = 0.079, v, z = 0.0395, t0 = 0.26, sv = 0.063, s = 0.1)
model <- NH_speed_model  %>% group_by(response) %>% 
  summarise(prop = mean(response == "upper"), mean_rt = 1000*(mean(rt)), 
            median_rt = mean(rt), sd_rt = 2*(sd(rt))) %>% 
  ungroup()

xyplot(mean_rt ~prop, model, group = response, type = "b", 
       auto.key = list(lines = FALSE), col = c("darkorchid1", "deepskyblue2"),
       xlab = "Response Probability",ylab = "Mean Reaction Time (ms)")

# plotting quantiles for NH speed condition
quantiles_data_NH <- NH_speed_model %>%
  group_by(response) %>%
  nest() %>% 
  mutate(quantiles = map(data, ~ as.data.frame(t(quantile(.x$rt, probs = quantiles))))) %>% 
  unnest(quantiles) %>% 
  gather("quantile_NH_speed", "rt",`10%`:`90%`) %>% 
  arrange(response)

NH_speed_quantile_model <- quantiles_data_NH  %>% group_by(quantile_NH_speed,response) %>% 
 summarise(prop = mean(response == "upper")
            , mean_rt = 1000*(mean(rt)
           ) 
   , median_rt = mean(rt)
  ) %>% 
 ungroup()

xyplot(mean_rt ~ prop, NH_speed_quantile_model, group = quantile_NH_speed, 
       type = "b", 
       auto.key = list(lines = FALSE, space = "right"), ylab = "Mean Reaction Time (ms)")

NH = data.frame(0)
for (i in 1:33) {
  NH[i] <- rdiffusion(100, a = 0.079, v = 0.9*(probability[i] - 0.45), z = 0.0395, t0 = 0.26, sv = 0.063, s = 0.1)
  
  print(NH[i])
}
# upper / lower = correct and incorrect responses
quantiles_NH_speed <- NH_speed_model %>%
  group_by(response) %>%
  nest() %>% 
  mutate(quantiles = map(data, ~ as.data.frame(t(quantile(.x$rt, probs = quantiles))))) %>% 
  unnest(quantiles) %>% 
  gather("quantile", "rt",`10%`:`90%`) %>% 
  arrange(response)

xyplot(rt*1000 ~ response, quantiles_NH_speed, group = quantile, type = "b", 
       auto.key = list(lines = FALSE, space = "right"), ylab = "Mean Reaction Time (ms)",ylim = c(200,600))


#data_1 <- data  %>% group_by(strength,id,instruction) %>% 
 # summarise(prop = mean(response == "dark"), mean_rt = 1000*(mean(rt)), 
  #          median_rt = mean(rt), sd_rt = 2*(sd(rt))) %>% 
  #ungroup()

#xyplot(mean_rt ~prop|id, data_1, group = instruction, type = "b", 
 #      auto.key = list(lines = FALSE), col = c("darkorchid1", "deepskyblue2"), 
  #     xlab = "Response Probability",ylab = "Mean Reaction Time (ms)",  
#)

# quantiles for NH accuracy model
NH_accuracy_model <- rdiffusion(4187, a = 0.16, v, z = 0.08, t0 = 0.26, sv = 0.063, s = 0.1)

quantiles_NH_accuracy <- NH_accuracy_model %>%
  group_by(response) %>%
  nest() %>% 
  mutate(quantiles = map(data, ~ as.data.frame(t(quantile(.x$rt, probs = quantiles))))) %>% 
  unnest(quantiles) %>% 
  gather("quantile", "rt",`10%`:`90%`) %>% 
  arrange(response)

xyplot(rt*1000 ~ response, quantiles_NH_accuracy, group = quantile, type = "b", 
       auto.key = list(lines = FALSE, space = "right"), ylab = "Mean Reaction Time (ms)",ylim = c(100,2000))

#NH = data.frame(quantiles_NH_speed$rt, q_diff)
#plot(NH)
###### Histogram of RT data (Fig. 5)
NH_RT = data$rt[data$instruction == "accuracy"]
RT = NH_RT[data$rt]
RT_acc = data[data$rt & data$instruction == "accuracy"]
NH_RT_acc_hist = hist(data$rt[data$instruction == "accuracy" & data$id == "nh" & data$correct == TRUE], breaks = 60, xlim = c(0,1.5), xlab = "Time (s)")
NH_RT_speed_hist = hist(data$rt[data$instruction == "speed" & data$id == "nh" & data$correct == TRUE], breaks = 50, xlim = c(0,1.5), xlab = "Time (s)")

######## JF #######
n_JF_acc = rr98[rr98$id == "jf" & rr98$instruction == "accuracy",] # number of accuracy trials
n_JF_speed = rr98[rr98$id == "jf" & rr98$instruction == "speed",] # number of speed trials

JF_speed_model <- rdiffusion(3945, a = 0.08, v, z = 0.04, t0 = 0.274, sv = 0.093, s = 0.1)
JF_accuracy_model <- rdiffusion(3943, a = 0.148, v, z = 0.074, t0 = 0.274, sv = 0.093, s = 0.1)

## Speed Model
quantiles_JF_speed <- JF_speed_model %>%
  group_by(response) %>%
  nest() %>% 
  mutate(quantiles = map(data, ~ as.data.frame(t(quantile(.x$rt, probs = quantiles))))) %>% 
  unnest(quantiles) %>% 
  gather("quantile", "rt",`10%`:`90%`) %>% 
  arrange(response)

xyplot(rt*1000 ~ response, quantiles_JF_speed, group = quantile, type = "b", 
       auto.key = list(lines = FALSE, space = "right"), ylab = "Mean Reaction Time (ms)",ylim = c(200,600))

## Accuracy Model

quantiles_JF_accuracy <- JF_accuracy_model %>%
  group_by(response) %>%
  nest() %>% 
  mutate(quantiles = map(data, ~ as.data.frame(t(quantile(.x$rt, probs = quantiles))))) %>% 
  unnest(quantiles) %>% 
  gather("quantile", "rt",`10%`:`90%`) %>% 
  arrange(response)

xyplot(rt*1000 ~ response, quantiles_JF_accuracy, group = quantile, type = "b", 
       auto.key = list(lines = FALSE, space = "right"), ylab = "Mean Reaction Time (ms)",ylim = c(100,2000))

######## KR ##########
n_KR_acc = rr98[rr98$id == "kr" & rr98$instruction == "accuracy",] # number of accuracy trials
n_KR_speed = rr98[rr98$id == "kr" & rr98$instruction == "speed",] # number of speed trials

KR_speed_model <- rdiffusion(3859, a = 0.073, v, z = 0.0365, t0 = 0.228, sv = 0.082, s = 0.1)
KR_accuracy_model <- rdiffusion(3921, a = 0.186, v, z = 0.093, t0 = 0.228, sv = 0.082, s = 0.1)

## Speed Model
quantiles_KR_speed <- KR_speed_model %>%
  group_by(response) %>%
  nest() %>% 
  mutate(quantiles = map(data, ~ as.data.frame(t(quantile(.x$rt, probs = quantiles))))) %>% 
  unnest(quantiles) %>% 
  gather("quantile", "rt",`10%`:`90%`) %>% 
  arrange(response)

xyplot(rt*1000 ~ response, quantiles_KR_speed, group = quantile, type = "b", 
       auto.key = list(lines = FALSE, space = "right"), ylab = "Mean Reaction Time (ms)",ylim = c(200,600))

## Accuracy Model
quantiles_KR_accuracy <- KR_accuracy_model %>%
  group_by(response) %>%
  nest() %>% 
  mutate(quantiles = map(data, ~ as.data.frame(t(quantile(.x$rt, probs = quantiles))))) %>% 
  unnest(quantiles) %>% 
  gather("quantile", "rt",`10%`:`90%`) %>% 
  arrange(response)

xyplot(rt*1000 ~ response, quantiles_KR_accuracy, group = quantile, type = "b", 
       auto.key = list(lines = FALSE, space = "right"), ylab = "Mean Reaction Time (ms)", ylim = c(0,2000))
