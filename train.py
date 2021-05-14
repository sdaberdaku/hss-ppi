##################################################
#### The following code is used to ensure the ####
#### reproducibility of the results           ####
##################################################
# Setting a seed value to ensure result reproducibility
seed_value = 43

import os
os.environ['PYTHONHASHSEED']=str(seed_value)

import random
random.seed(seed_value)

import numpy as np
np.random.seed(seed_value)

import tensorflow as tf
tf.compat.v1.set_random_seed(seed_value)
##################################################

from tensorflow.keras import optimizers

from sklearn.metrics import roc_curve, auc, precision_recall_curve, f1_score, accuracy_score, precision_score, recall_score, matthews_corrcoef

import matplotlib.pyplot as plt
import pandas as pd

from utils.utils import *
from utils.models import GCN, MLP
from utils.config import args
from utils.load_protein_data import load_protein_data

import gzip
import pickle
from pathlib import Path

# apply threshold to positive probabilities to create labels
def to_labels(current_preds, threshold):
	return (current_preds >= threshold).astype('int')


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
print('tf version:', tf.__version__)
assert tf.__version__.startswith('2.')

import warnings
warnings.filterwarnings("ignore")

# load data

protein_class = args.protein_class
is_bound = args.is_bound
is_receptor = args.is_receptor
distance_threshold = args.distance_threshold
representation_type = args.representation_type


lr = "r" if is_receptor else "l"
bu = "b" if is_bound else "u"
dataset_pickle_filename = f"data/{protein_class}_{lr}_{bu}.pkl.gz"
if args.load_pickle :
    with gzip.open(dataset_pickle_filename, 'rb') as file_handler:
        adj, features, y_train, y_val, y_test, train_mask, val_mask, test_mask, train_groups, val_groups, test_groups = pickle.load(file_handler)
else :
    adj, features, y_train, y_val, y_test, train_mask, val_mask, test_mask, train_groups, val_groups, test_groups = load_protein_data(aa_index_file="data/AAIndexFeatures.csv", representation_type=representation_type, protein_class=protein_class, is_bound=is_bound, is_receptor=is_receptor, distance_threshold=distance_threshold)
    with gzip.open(dataset_pickle_filename, 'wb') as file_handler:
        pickle.dump([adj, features, y_train, y_val, y_test, train_mask, val_mask, test_mask, train_groups, val_groups, test_groups], file_handler)

print('adj:', adj.shape)
print('features:', features.shape)
print('y:', y_train.shape, y_val.shape, y_test.shape)
print('mask:', train_mask.shape, val_mask.shape, test_mask.shape)

# D^-1@X
features = preprocess_features(features)  # [49216, 2], [49216], [2708, 1433]

print('features coordinates::', features[0].shape)
print('features data::', features[1].shape)
print('features shape::', features[2])

support = [preprocess_adj(adj)]
num_supports = 1
model_func = GCN

# Create model
model = GCN(input_dim=features[2][1],
            output_dim=y_train.shape[1],
            num_features_nonzero=features[1].shape)

train_label = tf.convert_to_tensor(y_train)
train_mask = tf.convert_to_tensor(train_mask)
val_label = tf.convert_to_tensor(y_val)
val_mask = tf.convert_to_tensor(val_mask)
test_label = tf.convert_to_tensor(y_test)
test_mask = tf.convert_to_tensor(test_mask)
features = tf.SparseTensor(*features)

support = [tf.cast(tf.SparseTensor(*support[0]), dtype=tf.float32)]
num_features_nonzero = features.values.shape
dropout = args.dropout

#optimizer = optimizers.RMSprop(learning_rate=args.learning_rate)
optimizer = optimizers.Adam(learning_rate=args.learning_rate)

results_path = f"results(epochs={args.epochs},representation_type={representation_type},distance_threshold={distance_threshold})"        
Path(results_path).mkdir(parents=True, exist_ok=True)

output_filename = f"{results_path}/{protein_class}_{lr}_{bu}_interface_prediction.txt"
with open(output_filename, "w") as out_file:
    out_file.write("Epoch\tloss_train\tAUC_train\tloss_val\tAUC_val\tloss_test\tAUC_test\n")

for epoch in range(args.epochs):

    with tf.GradientTape() as tape:
        loss_train, auc_train, _ = model((features, train_label, train_mask, support))
        grads = tape.gradient(loss_train, model.trainable_variables)
        optimizer.apply_gradients(zip(grads, model.trainable_variables))

    loss_val, auc_val, val_pred = model((features, val_label, val_mask, support), training=False)
    loss_test, auc_test, _ = model((features, test_label, test_mask, support), training=False)

    print("Epoch:\t", epoch,
          "loss_train:\t", float(loss_train),
          "AUC_train:\t", float(auc_train),
          "loss_val:\t", float(loss_val),
          "AUC_val:\t", float(auc_val),
          "loss_test:\t", float(loss_test),
          "AUC_test:\n", float(auc_test))

    with open(output_filename, "a") as out_file:
        out_file.write("%d\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n" %
                        (epoch, loss_train, auc_train, loss_val, auc_val, loss_test, auc_test))

print("\n\nComputing the optimal classification threshold on the validation set.\n\n")
val_pdb_ids = sorted(list(set(val_groups)))
    
val_lab = tf.boolean_mask(val_label, val_mask).numpy().flatten()

step=0.01 # classification threshold step
thresholds = np.arange(0, 1, step)
thresholds_matrix = np.empty((0, len(thresholds)))
print("Validation set: ", val_pdb_ids)

for pdb_id in val_pdb_ids:
    indices = [i for i, x in enumerate(val_groups) if x == pdb_id]
    current_preds_val = [val_pred[i] for i in indices]
    current_labels_val = [val_lab[i] for i in indices]
    
    # evaluate each threshold
    current_thresholds= np.array([f1_score(current_labels_val, to_labels(current_preds_val, t)) for t in thresholds])  
    thresholds_matrix = np.vstack((thresholds_matrix, current_thresholds))

    
print(thresholds_matrix)
average_thresholds=np.mean(thresholds_matrix,axis=0)

print(average_thresholds)

selected_threshold_index = np.argmax(average_thresholds)      
print("Selected threshold index: ", selected_threshold_index)
selected_threshold = (selected_threshold_index + 1) * step
print("Selected threshold: ", selected_threshold)

test_loss, test_auc, test_pred = model((features, test_label, test_mask, support), training=False)

print("test_loss:", float(test_loss), "\ttest_auc:", float(test_auc))

test_pdb_ids = sorted(list(set(test_groups)))
print("Test set: ", test_pdb_ids)

fpr = dict()
tpr = dict()
roc_auc = dict()

y_lab = tf.boolean_mask(test_label, test_mask).numpy().flatten()

f, t, _ = roc_curve(y_lab, test_pred)
output_roc_filename = f"{results_path}/{protein_class}_{lr}_{bu}_roc_plot_mic.pdf"
plt.figure()
plt.step(f, t, where='post')
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.ylim([0.0, 1.0])
plt.xlim([0.0, 1.0])
plt.savefig(output_roc_filename)

pred_out_dir = f"{results_path}/predictions/test/{protein_class}"
if not os.path.exists(pred_out_dir):
    os.makedirs(pred_out_dir)

score_tot = 0
acc_tot = 0
prec_tot = 0
rec_tot = 0
mcc_tot = 0
for pdb_id in test_pdb_ids:
    indices = [i for i, x in enumerate(test_groups) if x == pdb_id]
    current_preds = [test_pred[i] for i in indices]
    current_labels = [y_lab[i] for i in indices]     

    # evaluate each threshold calcolato per f1 sulle molecole della validation_set    
    scores = f1_score(current_labels, to_labels(current_preds, selected_threshold))    
    acc = accuracy_score(current_labels, to_labels(current_preds, selected_threshold))
    prec = precision_score(current_labels, to_labels(current_preds, selected_threshold)) 
    rec = recall_score(current_labels, to_labels(current_preds, selected_threshold))
    mcc= matthews_corrcoef(current_labels, to_labels(current_preds, selected_threshold))
    
    
    with open(f"{pred_out_dir}/{pdb_id}_pred.txt", "w") as out_pred:
        out_pred.write("\n".join(str(i) for i in current_preds))
          
    fpr[pdb_id], tpr[pdb_id], _ = roc_curve(current_labels, current_preds)
    roc_auc[pdb_id] = auc(fpr[pdb_id], tpr[pdb_id])    
    print("###################")
    print('AU-ROC:%s %.5f' %(pdb_id,roc_auc[pdb_id]))
    print('Threshold=%.3f, F-Score=%.5f, acc=%.5f, prec=%.5f, rec=%.5f' % (selected_threshold, scores, acc, prec,rec))
    score_tot += scores
    acc_tot += acc
    prec_tot += prec
    rec_tot += rec
    mcc_tot += mcc
    
avg_roc_auc = sum(roc_auc.values()) / len(roc_auc)
score_tot = score_tot/len(roc_auc)
acc_tot = acc_tot/len(roc_auc)
prec_tot = prec_tot/len(roc_auc)
rec_tot = rec_tot/len(roc_auc)
mcc_tot = mcc_tot/len(roc_auc)
 
print("Macro-Average AU-ROC: %.5f" % avg_roc_auc)
print("Macro-Average score: %.5f" % score_tot)
print("Macro-Average acc_tot: %.5f" % acc_tot)
print("Macro-Average prec_tot: %.5f" % prec_tot)
print("Macro-Average rec_tot: %.5f" % rec_tot)
print("Macro-Average mcc_tot: %.5f" % mcc_tot)

p, r, _ = precision_recall_curve(y_true=y_lab, probas_pred=test_pred)
output_pr_filename = f"{results_path}/{protein_class}_{lr}_{bu}_pr_plot_mic.pdf"
plt.figure()
plt.step(r, p, where='post')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.0])
plt.xlim([0.0, 1.0])

plt.savefig(output_pr_filename)
