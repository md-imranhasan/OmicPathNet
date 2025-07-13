# === train_hub_ranker_cv.R ===
library(xgboost)
library(Matrix)
library(caret)
library(ggplot2)

set.seed(42)

# === Load your real features
df <- read.csv("ppi_hub_features.csv")

# Check
print(head(df))

# === Features and label ===
X <- as.matrix(df[, c("Degree", "Betweenness", "Closeness")])
y <- df$HubScore

# === Cross-validation folds ===
folds <- createFolds(y, k = 5, list = TRUE, returnTrain = FALSE)

dtrain <- xgb.DMatrix(data = X, label = y)

# === Params ===
params <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse",
  max_depth = 3,
  eta = 0.1
)

# === CV ===
cv <- xgb.cv(
  params = params,
  data = dtrain,
  nrounds = 100,
  nfold = 5,
  verbose = 1,
  early_stopping_rounds = 10,
  showsd = TRUE
)

print(cv)

best_nrounds <- cv$best_iteration
cat("✅ Best nrounds:", best_nrounds, "\n")

# === Final train ===
model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = best_nrounds
)

# === Save model ===
saveRDS(model, "hub_ranker_cv.rds")
cat("✅ Saved: hub_ranker_cv.rds\n")

# === Feature importance ===
imp <- xgb.importance(model = model)
print(imp)

# Plot importance
xgb.plot.importance(imp)

# Save importance as CSV
write.csv(imp, "hub_ranker_feature_importance.csv", row.names = FALSE)

# Optional: Save a pretty ggplot
imp_df <- as.data.frame(imp)
ggplot(imp_df, aes(x = reorder(Feature, Gain), y = Gain)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(x = "Feature", y = "Importance (Gain)", title = "XGBoost Feature Importance") +
  theme_minimal(base_size = 14)
ggsave("hub_ranker_feature_importance.png", width = 8, height = 5, dpi = 300)
