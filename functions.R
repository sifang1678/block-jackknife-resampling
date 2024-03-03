

mvmr_using_prs <- function(exposure_list, exposure_prs_list, outcome_list, covar_list, data) {
  l <- list()
  for (i in 1:length(outcome_list)) {
    out <- as.data.frame(matrix(data = NA, nrow = 4))
    outcome <- outcome_list[i]
    exposure1 <- exposure_list[1]
    exposure2 <- exposure_list[2]
    exposure_prs1 <- exposure_prs_list[1]
    exposure_prs2 <- exposure_prs_list[2]
    
    print(outcome)
    
    out[,1] <- outcome
    out[,2] <- rep(c(exposure1, exposure2), 2)
    out[,3] <- c("univariate", "univariate", "multivariate", "multivariate")
    
    ivmodel1 <- ivreg(data = data, 
                      as.formula(paste0("rntransform(`", outcome, "`) ~ ", exposure1, " + ", paste(covar_list, collapse = "+"), " | ", 
                                        exposure_prs1, " + ", paste(covar_list, collapse = "+"))))
    ivmodel2 <- ivreg(data = data, 
                      as.formula(paste0("rntransform(`", outcome, "`) ~ ", exposure2, " + ", paste(covar_list, collapse = "+"), " | ", 
                                        exposure_prs2, " + ", paste(covar_list, collapse = "+"))))
    ivmodel3 <- ivreg(data = data, 
                      as.formula(paste0("rntransform(`", outcome, "`) ~ ", exposure1, " + ", exposure2, " + ", paste(covar_list, collapse = "+"), " | ", 
                                        exposure_prs1, " + ", exposure_prs2, " + ", paste(covar_list, collapse = "+"))))
    res1 <- summary(ivmodel1, vcov = sandwich, df = Inf, diagnostics = TRUE)
    res2 <- summary(ivmodel2, vcov = sandwich, df = Inf, diagnostics = TRUE)
    res3 <- summary(ivmodel3, vcov = sandwich, df = Inf, diagnostics = TRUE)
    
    out[,c(4:7)] <- rbind(res1$coef[2,], res2$coef[2,]) %>% rbind(res3$coef[c(2,3),])
    out[,c(8:9)] <- rbind(res1$diagnostics[1, 3:4], res2$diagnostics[1, 3:4]) %>% rbind(res3$diagnostics[c(1,2), 3:4])
    
    names(out) <- c("outcome", "exposure", "model", "b", "se", "tval", "pval", "F", "F_pval")
    
    l[[i]] <- out
  }
  
  l_df <- do.call(rbind, l)
  return(l_df)
}


mvmr_using_prs_binary <- function(exposure_list, exposure_prs_list, outcome_list, covar_list, data) {
  l <- list()
  for (i in 1:length(outcome_list)) {
    out <- as.data.frame(matrix(data = NA, nrow = 4))
    outcome <- outcome_list[i]
    exposure1 <- exposure_list[1]
    exposure2 <- exposure_list[2]
    exposure_prs1 <- exposure_prs_list[1]
    exposure_prs2 <- exposure_prs_list[2]
    
    print(outcome)
    
    attach(data)
    
    con = data %>% dplyr::filter((!!as.name(outcome)) == 0) %>% dplyr::select(all_of(exposure_prs1)) %>% unlist #values for prs in the controls only
    lm.1 <- lm(as.formula(paste0(exposure1, "~", exposure_prs1, " + ", paste(covar_list, collapse = "+"))), 
               data = data %>% dplyr::filter((!!as.name(outcome)) == 0))
    predict.con.1 = predict(lm.1, newdata=list(con=data %>% dplyr::select(all_of(exposure_prs1)))) #Generate predicted values for all participants based on the linear regression in the controls only.  
    tsls.con.1 = glm(as.formula(paste0(outcome, "~ predict.con.1 + ", paste(covar_list, collapse = "+"))),
                     data = data, family="binomial") #Fit a logistic regression model on all the participants
    f.1 = summary(lm.1)$f[1]
    
    con = data %>% dplyr::filter((!!as.name(outcome)) == 0) %>% dplyr::select(all_of(exposure_prs2)) %>% unlist #values for prs in the controls only
    lm.2 <- lm(as.formula(paste0(exposure2, "~", exposure_prs2, " + ", paste(covar_list, collapse = "+"))), 
               data = data %>% dplyr::filter((!!as.name(outcome)) == 0))
    predict.con.2 = predict(lm.2, newdata=list(con=data %>% dplyr::select(all_of(exposure_prs2)))) #Generate predicted values for all participants based on the linear regression in the controls only.  
    tsls.con.2 = glm(as.formula(paste0(outcome, "~ predict.con.2 + ", paste(covar_list, collapse = "+"))),
                     data = data, family="binomial") #Fit a logistic regression model on all the participants
    f.2 = summary(lm.2)$f[1]
    
    tsls.con.mvmr = glm(as.formula(paste0(outcome, "~ predict.con.1 + predict.con.2 + ", paste(covar_list, collapse = "+"))),
                        data = data, family="binomial")
    
    est <- felm(as.formula(paste0(outcome, "~ ", paste(covar_list, collapse = "+"), "| 0 | (", exposure1, " | ", exposure2, " ~ ", exposure_prs1, " + ", exposure_prs2, ")")),
                data = data)
    f.34 <- condfstat(est)[1:2]
    
    res1 <- summary(tsls.con.1)
    res2 <- summary(tsls.con.2)
    res3 <- summary(tsls.con.mvmr)
    
    out[,1] <- outcome
    out[,2] <- rep(c(exposure1, exposure2), 2)
    out[,3] <- c("univariate", "univariate", "multivariate", "multivariate")
    out[,c(4:7)] <- rbind(res1$coef[2,], res2$coef[2,]) %>% rbind(res3$coef[c(2,3),])
    out[,8] <- c(f.1, f.2, f.34)
    
    names(out) <- c("outcome", "exposure", "model", "b", "se", "zval", "pval", "F")
    
    l[[i]] <- out
    
    detach(data)
  }
  
  l_df <- do.call(rbind, l)
  return(l_df)
}


mr_using_prs <- function(exposure, exposure_prs, outcome_list, covar_list, data) {
  l <- list()
  for (i in 1:length(outcome_list)) {
    out <- as.data.frame(matrix(data = NA, nrow = 1))
    outcome <- outcome_list[i]
    
    print(outcome)
    
    out[,1] <- outcome
    out[,2] <- exposure
    out[,3] <- "univariate"
    
    ivmodel <- ivreg(data = data, 
                     as.formula(paste0("`", outcome, "` ~ ", exposure, " + ", paste(covar_list, collapse = "+"), " | ", 
                                       exposure_prs, " + ", paste(covar_list, collapse = "+"))))
    res <- summary(ivmodel, vcov = sandwich, df = Inf, diagnostics = TRUE)
    
    out[,4:7] <- res$coef[2,]
    out[,c(8:9)] <- res$diagnostics[1, 3:4]
    
    names(out) <- c("outcome", "exposure", "model", "b", "se", "tval", "pval", "F", "F_pval")
    
    l[[i]] <- out
  }
  
  l_df <- do.call(rbind, l)
  return(l_df)
}


mr_using_prs_binary <- function(exposure, exposure_prs, outcome_list, covar_list, data) {
  l <- list()
  for (i in 1:length(outcome_list)) {
    out <- as.data.frame(matrix(data = NA, nrow = 1))
    outcome <- outcome_list[i]
    
    print(outcome)
    
    tmp <- data %>% dplyr::select(all_of(c(exposure, exposure_prs, outcome_list, covar_list)))
    tmp <- na.omit(tmp)
    attach(tmp)
    
    con = tmp %>% dplyr::filter((!!as.name(outcome)) == 0) %>% dplyr::select(all_of(exposure_prs)) %>% unlist #values for prs in the controls only
    lm <- lm(as.formula(paste0("`", exposure, "`~`", exposure_prs, "` + ", paste(covar_list, collapse = "+"))), 
             data = tmp %>% dplyr::filter((!!as.name(outcome)) == 0))
    predict.con = predict(lm, newdata=list(con=tmp %>% dplyr::select(all_of(exposure_prs)))) #Generate predicted values for all participants based on the linear regression in the controls only.  
    tsls.con = glm(as.formula(paste0("`", outcome, "`~ predict.con + ", paste(covar_list, collapse = "+"))),
                   data = tmp, family="binomial") #Fit a logistic regression model on all the participants
    f = summary(lm)$f[1]
    
    res <- summary(tsls.con)
    
    out[,1] <- outcome
    out[,2] <- exposure
    out[,3] <- c("univariate")
    out[,c(4:7)] <- res$coef[2,]
    out[,8] <- f
    
    names(out) <- c("outcome", "exposure", "model", "b", "se", "zval", "pval", "F")
    
    l[[i]] <- out
    
    detach(tmp)
  }
  
  l_df <- do.call(rbind, l)
  return(l_df)
}




plot_mvmr <- function(data, outcome_labels, exposure_labels) {
  p <- ggplot(data, aes(x = exposure, y = b, color = exposure)) +
    geom_point() +
    geom_errorbar(aes(ymin = b - 1.96 * se, ymax = b + 1.96 * se), width = 0.2) +
    facet_grid(outcome ~ factor(model, levels = c("univariate", "multivariate")), 
               switch = "both",
               labeller = labeller(outcome = outcome_labels)) +
    scale_color_discrete(breaks = names(exposure_labels), 
                         labels = unname(exposure_labels)) +
    geom_hline(yintercept = 0, color = "blue", linetype = 2) +
    labs(title = "One-Sample MR", x = "",
         y = "Beta Coefficients (95% Confidence Interval)") +
    theme_bw() +
    theme(strip.text.y.left = element_text(angle = 0, face="bold"),
          plot.title = element_text(size=16, face="bold", hjust = 0.5),
          plot.subtitle = element_text(size=10, face="bold", hjust = 0.5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(face="bold"),
          axis.title = element_text(size=12, face="bold"),
          legend.position="bottom") +
    coord_flip()
  
  return(p)
}


plot_sex_stratified_mvmr <- function(data, outcome_labels, exposure_labels) {
  p <- ggplot(data, aes(x = exposure, y = b, color = exposure)) +
    geom_point() +
    geom_errorbar(aes(ymin = b - 1.96 * se, ymax = b + 1.96 * se), width = 0.2) +
    facet_grid(outcome ~ sex + factor(model, levels = c("univariate", "multivariate")), 
               switch = "both",
               labeller = labeller(outcome = outcome_labels)) +
    scale_color_discrete(breaks = names(exposure_labels), 
                         labels = unname(exposure_labels)) +
    geom_hline(yintercept = 0, color = "blue", linetype = 2) +
    labs(title = "One-Sample MR (sex-stratified)", x = "",
         y = "Beta Coefficients (95% Confidence Interval)") +
    theme_bw() +
    theme(strip.text.y.left = element_text(angle = 0, face="bold"),
          plot.title = element_text(size=16, face="bold", hjust = 0.5),
          plot.subtitle = element_text(size=10, face="bold", hjust = 0.5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(face="bold"),
          axis.title = element_text(size=12, face="bold"),
          legend.position="bottom") +
    coord_flip()
  
  return(p)
}


plot_mvmr_binary <- function(data, outcome_labels, exposure_labels) {
  data <- data %>%
    mutate(or = exp(b),
           ci_higher = exp(b + 1.96 * se),
           ci_lower = exp(b - 1.96 * se))
  p <- ggplot(data, aes(x = exposure, y = or, color = exposure)) +
    geom_point() +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_higher), width = 0.2) +
    facet_grid(outcome ~ factor(model, levels = c("univariate", "multivariate")), 
               switch = "both",
               labeller = labeller(outcome = outcome_labels)) +
    scale_color_discrete(breaks = names(exposure_labels), 
                         labels = exposure_labels) +
    geom_hline(yintercept = 1, color = "blue", linetype = 2) +
    labs(title = "One-Sample MR", x = "",
         y = "Odds Ratio (95% Confidence Interval)") +
    theme_bw() +
    theme(strip.text.y.left = element_text(angle = 0, face="bold"),
          plot.title = element_text(size=16, face="bold", hjust = 0.5),
          plot.subtitle = element_text(size=10, face="bold", hjust = 0.5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(face="bold"),
          axis.title = element_text(size=12, face="bold"),
          legend.position="bottom") +
    coord_flip()
  
  return(p)
}


plot_sex_stratified_mvmr_binary <- function(data, outcome_labels, exposure_labels) {
  data <- data %>%
    mutate(or = exp(b),
           ci_higher = exp(b + 1.96 * se),
           ci_lower = exp(b - 1.96 * se))
  p <- ggplot(data, aes(x = exposure, y = or, color = exposure)) +
    geom_point() +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_higher), width = 0.2) +
    facet_grid(outcome ~ sex + factor(model, levels = c("univariate", "multivariate")), 
               switch = "both",
               labeller = labeller(outcome = outcome_labels)) +
    scale_color_discrete(breaks = names(exposure_labels), 
                         labels = unname(exposure_labels)) +
    geom_hline(yintercept = 1, color = "blue", linetype = 2) +
    labs(title = "One-Sample MR (sex-stratified)", x = "",
         y = "Odds Ratio (95% Confidence Interval)") +
    theme_bw() +
    theme(strip.text.y.left = element_text(angle = 0, face="bold"),
          plot.title = element_text(size=16, face="bold", hjust = 0.5),
          plot.subtitle = element_text(size=10, face="bold", hjust = 0.5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(face="bold"),
          axis.title = element_text(size=12, face="bold"),
          legend.position="bottom") +
    coord_flip()
  
  return(p)
}


