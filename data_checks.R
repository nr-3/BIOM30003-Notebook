## Checks

check_mat <- as.matrix(monkeyHSPC@assays@data@listData$proportion)
check1 <- as.data.frame(check_mat)
check2 <- check1

check1 <- check1 |>
  pivot_longer(cols = all_of(colnames(check1)), names_to = 'sample', values_to = 'proportion') |>
  group_by(sample)

# Checking if proportions are truly adding up to 1
check1 |>
  summarise(sum = sum(proportion)) |>
  arrange(sum)


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
