## Checks

check_mat <- as.matrix(monkeyHSPC@assays@data@listData$proportion)
check1 <- as.data.frame(check_mat)
check2 <- check1

check1$barcode <- row.names(check1)

check1 <- check1 |>
  pivot_longer(cols = all_of(colnames(check1)[-43]), names_to = 'sample', values_to = 'proportion') |>
  group_by(sample)

# Checking if proportions are truly adding up to 1
check1 |>
  summarise(sum = sum(proportion)) |>
  arrange(sum)

check1 |>
  ungroup() |>
  group_by(barcode) |>
  summarise(zero_eff = n()) |>
  arrange(desc(zero_eff))

# Check distribution of the proportions
check2 <- check2[,monkeyHSPC$sampleMetadata$Celltype == 'NK_CD56n_CD16p']
check2 |>
  pivot_longer(cols = all_of(colnames(check2)), names_to = 'sample', values_to = 'proportion') |>
  ggplot(aes(x = factor(0), y = asin(sqrt(proportion)))) +
  ggbeeswarm::geom_quasirandom()

check3 <- check2[,monkeyHSPC$sampleMetadata$Celltype != 'NK_CD56n_CD16p']
check3 |>
  pivot_longer(cols = all_of(colnames(check3)), names_to = 'sample', values_to = 'proportion') |>
  ggplot(aes(x = factor(0), y = asin(sqrt(proportion)))) +
  ggbeeswarm::geom_quasirandom()

View(monkeyHSPC@assays@data@listData$proportion)
View(se_HSPC@assays)

# Average zero percentage

check1 |>
  summarise(prop = sum(proportion==0)/n()) |>
  summarise(mean = mean(prop))

# Checking implantation technique

s <- as.vector(simulation_se@assays@data@listData[[1]][100,])
original <- s
s <- rescale_beta(s[grp1],beta_diff,zero_effect)
simulation_se@assays@data@listData[[1]][100,grp1] <- s$rescaled
# Should be around 15
sum(simulation_se@assays@data@listData[[1]][100,] == original)


sim1_res |>
  group_by(model,beta) |>
  summarise(
    sd = sqrt(var(AUC)) / sqrt(n())
  )

sim_res |>
  group_by(model, beta) |>
  summarise(
    mean = mean(AUC),
    se = sqrt(var(AUC)) / sqrt(n()),
    .groups = "drop"
  )

sim1_res |>
  filter(model == "LM", beta == 0.7) |>
  pull(AUC) |>
  sd()


sim_res |>
  group_by(model, beta) |>
  summarise(AUC = mean(AUC), se = sd(AUC) / sqrt(n()))
