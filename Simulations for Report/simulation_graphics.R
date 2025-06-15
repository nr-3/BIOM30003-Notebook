######################## Simulation via implantation ################################

# Register parallel backend
num_cores <- parallel::detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

zero_effect <- 0
beta_params <- 0.95
num_sim <- 3
n_rows <- dim(se_HSPC_2grp)[[1]]
n_cols <- dim(se_HSPC_2grp)[[2]]
TP <- n_rows*0.1

#sim_res <- data.frame(matrix(ncol = length(colnames), nrow = 0))
#colnames(sim_res) <- colnames

for (beta_diff in beta_params) {
  res <- foreach(
    j = 1:num_sim,
    .combine = rbind,
    .export = c("all_sims", "rescale_beta", "lm_sim", "limma_sim", 
                "wilcoxon_t_sim", "strainspy_sim", "se_HSPC_2grp", 
                "n_cols", "n_rows", "colnames", "zero_effect", "beta_diff")
  ) %dopar% {
    all_sims(j, beta_diff, zero_effect)
  }
  
  sim_res <- rbind(sim_res, res)
  write.csv(sim_res, file = 'C:/Users/nraja/Downloads/sim2_res.csv', row.names = F)
}

## Adding Specificity
sim_res$specificity <- (n_rows - (TP + TP * (1 - sim_res$precision) / sim_res$precision + TP * (1 - sim_res$recall) / sim_res$recall)) /
  (n_rows - (TP + TP * (1 - sim_res$recall) / sim_res$recall))

stopCluster(cl)

t################################ Graphics ######################################
#sim_res[,c(4,5,6)] <- NA
#combined <- #rbind(sim_res_adj, sim_res)

# Simulation Results: 1 Predictor
p1 <- plot_AUC(sim1_res)
p2 <- plot_FDR(sim1_res)
p3 <- plot_recall(sim1_res)

onepred <- ggpubr::ggarrange(p1,p2,p3, labels = c('A','B','C'),
                  ncol = 2, nrow = 2, common.legend = T, legend = 'right')
#onepred <- ggpubr::annotate_figure(onepred, bottom = 'Percentage Modifier of Implant')

# Simulation Results: 2 Predictors
p1 <- plot_AUC(sim2_res)
p2 <- plot_FDR(sim2_res)
p3 <- plot_recall(sim2_res)
twopred <- ggpubr::ggarrange(p1,p2,p3, labels = c('A','B','C'), 
                             ncol = 2, nrow = 2, common.legend = T, legend = 'right')
#twopred <- ggpubr::annotate_figure(twopred, bottom = 'Percentage Modifier of Implant')

onepred
twopred

# Simulation Results table: 1 Predictor
results <- sim1_res |>
  group_by(model) |>
  summarise(AUC = mean(AUC),
            Precision = mean(precision, na.rm = T),
            Recall = mean(recall, na.rm = T)) |>
  arrange(desc(AUC))

colnames(results)[1] <- 'Model'

rempsyc::nice_table(results,
                    title = 'Table 1: Average Metrics Across All One-Predictor Simulations')

# Simulation Results table: 2 predictors
results <- sim2_res |>
  group_by(model) |>
  summarise(AUC = mean(AUC),
            Precision = mean(precision, na.rm = T),
            Recall = mean(recall, na.rm = T)) |>
  arrange(desc(AUC))

colnames(results)[1] <- 'Model'

rempsyc::nice_table(results,
                    title = 'Table 2: Average Metrics Across All Two-Predictor Simulations')
