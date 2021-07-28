import pandas as pd
import sklearn
from sklearn.ensemble import RandomForestClassifier

# Load in data
metrics = pd.read_csv("/mnt/AchTeraD/Documents/Projects/scCUTseq/data/cell_classifier/metrics_LH.csv")

# Train classifier
clf = RandomForestClassifier(n_estimators=100, max_depth=None, min_samples_split=2, n_jobs=4,random_state=None,verbose=0)
(correct,predictions,predicted_prob)=train_ml(clf,X,Y,cv,verbose=False,plot=False)
(fpr,tpr,thres_roc)=roc_curve(correct,predicted_prob)   
plt.plot(fpr,tpr,label='Randomforest')