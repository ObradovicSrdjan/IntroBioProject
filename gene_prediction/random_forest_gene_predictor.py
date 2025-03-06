import pandas as pd
from pandas import DataFrame
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics import accuracy_score

genes_data_frame: DataFrame = pd.read_csv("data/generated_data/genes_training_set.csv")
sequences = genes_data_frame["Sequence"].values
labels = genes_data_frame["IsGene"].values

k = 3


def get_kmers(seq, k):
    return [seq[i : i + k] for i in range(len(seq) - k + 1)]


kmers_list = [" ".join(get_kmers(seq, k)) for seq in sequences]

vectorizer = CountVectorizer()
X = vectorizer.fit_transform(kmers_list)

X_train, X_test, y_train, y_test = train_test_split(
    X, labels, test_size=0.15, random_state=42
)

random_forest_classifier = RandomForestClassifier(n_estimators=200, random_state=42)
random_forest_classifier.fit(X_train, y_train)

y_pred = random_forest_classifier.predict(X_test)
print("Accuracy:", accuracy_score(y_test, y_pred))
