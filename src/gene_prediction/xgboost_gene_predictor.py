import pandas as pd
from pandas import DataFrame
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics import accuracy_score
from sklearn.model_selection import GridSearchCV, train_test_split
from xgboost import XGBClassifier

# Load the data
genes_data_frame: DataFrame = pd.read_csv("data/generated_data/genes_training_set.csv")
sequences = genes_data_frame["Sequence"].values
labels = genes_data_frame["IsGene"].values

k = 3


def get_kmers(seq, k):
    return [seq[i : i + k] for i in range(len(seq) - k + 1)]


# Convert sequences to k-mers
kmers_list = [" ".join(get_kmers(seq, k)) for seq in sequences]

# Convert k-mers to a numerical feature matrix
vectorizer = CountVectorizer()
X = vectorizer.fit_transform(kmers_list)

# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(
    X, labels, test_size=0.15, random_state=42
)

# Define the XGBoost classifier
xgb_classifier = XGBClassifier(random_state=42)

# Define the parameter grid
param_grid = {
    "n_estimators": [100, 200, 300],
    "max_depth": [3, 5, 7],
    "learning_rate": [0.01, 0.1, 0.2],
    "subsample": [0.8, 1.0],
    "colsample_bytree": [0.8, 1.0],
}

# Create the GridSearchCV object
grid_search = GridSearchCV(
    estimator=xgb_classifier, param_grid=param_grid, cv=3, n_jobs=-1, verbose=2
)

# Perform grid search
grid_search.fit(X_train, y_train)

# Get the best parameters and best estimator
best_params = grid_search.best_params_
best_estimator = grid_search.best_estimator_

print("Best Parameters:", best_params)

# Make predictions with the best estimator
y_pred = best_estimator.predict(X_test)

# Evaluate accuracy
print("Accuracy:", accuracy_score(y_test, y_pred))
