ss <- read.csv("starvation_Etoh.csv")
ss$replicate <- as.character(ss$replicate)
ss$geno <- factor(ss$geno,levels = c("WT", "daf-2", "daf-16", "daf-16; daf-2"))

survival_plot<-ggplot()+
  geom_point(data=ss, aes(x = day, y = prop_alive, color = geno))+
  geom_smooth(data=ss, aes(x = day, y = prop_alive, color = geno), method = "glm", method.args = list(family = "quasibinomial"), level = 0)+
  theme_bw()+
   ggtitle("survival rate ~ time")+labs(x="Days of L1 arrest",y="Proportion alive") 
survival_plot

result <- unique(subset(ss, select = "geno"))
for (i in 1:length(result$geno)) {
  model <- glm(data = ss, subset = (ss$geno == result$geno[i]),
               formula = prop_alive ~ day,
               family = "quasibinomial")
  half_life = - model$coefficients[1]/model$coefficients[2]
  goodness_of_fit <- 1 - model$deviance/model$null.deviance
  result$half_life[i] <- half_life
  result$goodness_of_fit[i] <- goodness_of_fit
}  
result

#round half_life to 2 digits
result$half_life <- round(result$half_life, 2)

result2 <- unique(subset(ss, select = c("geno", "replicate")))
result2$condition <- paste(result2$geno, result2$rep, sep = "_")

#then generate model
for (i in 1:length(result2$condition)) {
  model <- glm(data = ss, subset = (ss$geno == result2$geno[i]
                                    & ss$rep == result2$rep[i]),
               formula = prop_alive ~ day,
               family = "quasibinomial")
  half_life = - model$coefficients[1]/model$coefficients[2]
  goodness_of_fit <- 1 - model$deviance/model$null.deviance
  result2$half_life[i] <- half_life
  result2$goodness_of_fit[i] <- goodness_of_fit
}  
result2

result2$half_life <- round(result2$half_life, 2)
result2$rep <- as.factor(result2$rep)
bartlett.test(half_life~geno, data = result2)
pairwise.t.test(result2$half_life,result2$geno, pool.sd = TRUE, p.adjust = "none")


daf2xdaf16 <- result2
daf2xdaf16$geno %<>% as.character() %>% as.factor()
daf2xdaf16$strain %>% levels()

daf2xdaf16$daf16 <- NA
daf2xdaf16$daf2 <- NA
library("stringr")
daf2xdaf16[str_detect(daf2xdaf16$geno, "daf-16"), "daf16"] <- "Mut"
daf2xdaf16[str_detect(daf2xdaf16$geno, "daf-2"), "daf2"] <- "Mut"
daf2xdaf16$daf16[is.na(daf2xdaf16$daf16)] <- "WT"
daf2xdaf16$daf2[is.na(daf2xdaf16$daf2)] <- "WT"

# test for homogeneity of variance across groups
leveneTest(data = daf2xdaf16, half_life ~ daf16*daf2) # Pr 0.416, equal variance, can do ANOVA

bartlett.test(data = daf2xdaf16, half_life ~ interaction(daf16, daf2)) # Pr 0.0996, equal variance, can do ANOVA

# ANOVA assumes normal distribution and equal variance
aov(data = daf2xdaf16, half_life ~ daf16*daf2) %>% summary() # daf2:daf16 interaction 0.016
