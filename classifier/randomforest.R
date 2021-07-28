## Author: Luuk Harbers
## Date: xxxx-xx-xx
## Script for 

## Load/install packages
packages = c("data.table", "randomForest", "pROC")
sapply(packages, require, character.only = T)
source("/mnt/AchTeraD/Documents/R-functions/save_and_plot.R")

# Load in annotated metrics
metrics = fread("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/cell_classifier/metrics_LH.csv")

# Change character outcome to factor and only keep good/bad and remove sample/library columns
metrics[, user_quality := factor(user_quality, levels = c("good", "bad"))]
metrics = metrics[!is.na(user_quality), !c("sample", "library")]

size_training = 0.8

training_set = metrics[sample(.N, round(nrow(metrics) * size_training))]
validation_set = fsetdiff(metrics, training_set)

# Train model
rf = randomForest(user_quality ~ ., 
                  data = training_set, 
                  importance = T)

# Predict on validation set
prediction = as.data.table(predict(rf, validation_set, type="prob"))
prediction[, response := ifelse(good > bad, "good", "bad")]
prediction[, observation := validation_set$user_quality]
roc_curve = roc(prediction$observation, prediction$good)

roc_plt = ggroc(roc_curve, size = 1.2) +
  annotate("text", y = 0.1, x = 0.1, label = paste0("AUC = ", round(roc_curve$auc[[1]], 3)))

save_and_plot(roc_plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/classifier/roc_plot",
              width = 7, height = 7)

# Save randomforest
saveRDS(rf, "/mnt/AchTeraD/Documents/Projects/scCUTseq/data/cell_classifier/randomforest.rds")

# Get variable iomportance values
var_imps = data.frame(rf$importance)
var_imps$feature = rownames(var_imps)
setDT(var_imps)
setorder(var_imps, MeanDecreaseAccuracy)
var_imps[, feature := factor(feature, levels = feature)]

# Plot variable importances
imp_plt = ggplot(var_imps, aes(x = MeanDecreaseAccuracy, y = feature)) +
  geom_segment(aes(xend = 0, yend = feature), size = 1.2) +
  geom_point(size = 4, color = "orange") +
  labs(y = "Feature", x = "Mean decrease in accuracy")

save_and_plot(imp_plt, "/mnt/AchTeraD/Documents/Projects/scCUTseq/Plots/classifier/feature_importance",
              width = 6, height = 8)

# ROCR measurements
pred = prediction(prediction$good, prediction$observation)
f_measure = performance(pred, "f")
f_measure = max(f_measure@y.values[[1]])

