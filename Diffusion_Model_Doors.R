# Applying Diffusion Model to 3 Participants in the Doors Task

library("rtdists")
library(dplyr)

data_doors = read.csv("Doors_3.csv")
head(data_doors) # 1 = win, 2 = lose

data_doors_1 <- data_doors  %>% group_by(Participant, Trial) %>% 
  summarise(prop = mean(win_loss == 1), mean_rt = mean(RT), sd_rt = sd(RT)) %>% 
  ungroup()

xyplot(prop ~ Trial|Participant, data_doors_1, 
       type = "b", 
       auto.key = list(lines = TRUE), xlab = "Trial",
       ylab = "Proportion of Correct Responses")

# averages across participants
data_doors_2 <- data_doors  %>% group_by(Trial) %>% 
  summarise(prop = mean(win_loss == 1), mean_rt = mean(RT), sd_rt = sd(RT)) %>% 
  ungroup()

xyplot(prop ~ Trial, data_doors_2, 
       #group = instruction, 
       type = "b", 
       auto.key = list(lines = TRUE), xlab = "Trial",
       ylab = "Proportion of Correct Responses")

# Dividing the data into quantiles
quantiles <- c(0.1, 0.3, 0.5, 0.7, 0.9)
quantiles_data_doors <- data_doors %>%
  group_by(Participant, Trial) %>%
  nest() %>% 
  mutate(quantiles = map(data, ~ as.data.frame(t(quantile(.x$RT, probs = quantiles))))) %>% 
  unnest(quantiles) %>% 
  gather("quantile", "rt",`10%`:`90%`) %>% 
  arrange(Participant, Trial)

xyplot(rt ~ Trial|Participant, quantiles_data_doors, group = quantile, type = "b", 
       auto.key = list(lines = FALSE, space = "right"), ylab = "Mean Reaction Time (ms)")

data_doors_fitted <- data_doors %>% 
  group_by(Participant, Trial) %>% # loops across participant and trial
  nest()
data_doors_fitted

# function to assign drift rates to strengths 
objective_diffusion_separate <- function(pars, RT, win_los, drift, ...) {
  non_v_pars <- grep("^v", names(pars), invert = TRUE, value = TRUE)
  base_par <- length(non_v_pars)  # number of non-drift parameters
  densities <- vector("numeric", length(rt))
  for (i in seq_along(levels(drift))) {
    densities[drift == levels(drift)[i]] <- 
      ddiffusion(rt[drift == levels(drift)[i]], response=win_loss[drift == levels(drift)[i]], 
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
                                       rt = data$RT, response = data$win_loss, 
                                       drift = factor(data$Trial, seq_len(n_drift)))
      }
      cat("\nstart fitting.\n") 
      
      fit <- nlminb(start, objective_function, 
                    rt = data$RT, response = data$win_loss, 
                    drift = factor(data, seq_len(n_drift)))
      
      if (fit$objective < best_fit$objective) best_fit <- fit
    }
    out <- as.data.frame(t(unlist(best_fit[1:3])))
    colnames(out) <- sub("par.", "", colnames(out))
    out
  }

fit_diffusion <- data_doors_fitted %>% 
  mutate(fit = 
           map(data, 
               ~ensure_fit(data = ., start_function = get_start, 
                           objective_function = objective_diffusion_separate, 
                           as.numeric(response) %in% 1:2!!,
                           base_pars = c("a", "t0", "sv", "sz", "z")))) %>% 
  unnest(fit)




