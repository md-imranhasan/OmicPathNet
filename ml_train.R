# ml_train.R

library(xgboost)
library(data.table)

# ======= Simulate example training data =======
set.seed(42)
n <- 100  # 100 hubs

train_df <- data.frame(
  Degree = sample(1:200, n, replace = TRUE),
  Betweenness = runif(n),
  DrugHits = sample(0:10, n, replace = TRUE),
  Label = sample(0:1, n, replace = TRUE)  # Binary hub importance
)

print(head(train_df))

# ======= Prepare DMatrix =======
dmatrix <- xgb.DMatrix(
  data = as.matrix(train_df[, c("Degree", "Betweenness", "DrugHits")]),
  label = train_df$Label
)

# ======= Train simple XGBoost model =======
params <- list(
  objective = "binary:logistic",
  eval_metric = "logloss"
)

model <- xgb.train(
  params = params,
  data = dmatrix,
  nrounds = 30
)

# ======= Save the model =======
saveRDS(model, "ml_hub_ranker_model.rds")
print("✅ Model saved: ml_hub_ranker_model.rds")

# ======= (Optional) View feature importance =======
importance <- xgb.importance(model = model)
print(importance)
xgb.plot.importance(importance)
