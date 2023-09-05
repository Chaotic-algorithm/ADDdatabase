import pandas as pd
import ast

# Disclaimer:
# The code provided in this repository is a **sample** for benchmarking A-domain specificity prediction algorithms.
# Here it contains only benchmarking of the SANDPUMA tool using a dataset of all A-domains.
# A similar benchmarking procedure was applied for other tools and both bacterial and fungal datasets.
# Please note that this code is still in development to enhance its usability for a broader user base.
# Updates are coming soon!

def perf_measure(ground_truth, prediction):
    """
    Calculate true positives (TP), true negative (TN), false positive (FP) and false negative (FN)
    :param ground_truth: a binary list of the ground truth
    :param prediction: a binary list of predicted values
    :return: return a tuple of calculated counts
    """
    TP = 0
    FP = 0
    TN = 0
    FN = 0

    for i in range(len(prediction)):
        if ground_truth[i] == prediction[i] == 1:
            TP += 1
        if prediction[i] == 1 and ground_truth[i] != prediction[i]:
            FP += 1
        if ground_truth[i] == prediction[i] == 0:
            TN += 1
        if prediction[i] == 0 and ground_truth[i] != prediction[i]:
            FN += 1

    return (TP, FP, TN, FN)


def replace_characters(string):
    """
    Modify substrate's name into a uniform string
    :param string: substrate's name
    :return: lowercase string without additional characters apart from numbers and letters
    """
    string = string.replace(' ', '')
    string = ''.join(e for e in string if e.isalnum() or e.isspace())
    string = string.lower()
    return string


def calculate_metrics(ground_truth, prediction):
    """
    Calculate evaluation metrics from the confusion matrix
    :param ground_truth: list of ground truth substrates
    :param prediction: list of predicted substrates
    :return: micro-averaged evaluation metrics
    """
    per_class = {}
    metrics_dict = {}
    # list of all unique amino acids that appear in the ground truth list
    amino_acids = set(ground_truth)
    aa_list = []

    for amino_acid in amino_acids:
        # Turn ground truth and predicted substrate into binary labels for each amino acid
        as_true_binary = [1 if label == amino_acid else 0 for label in ground_truth]
        y_pred_binary = [1 if amino_acid in label else 0 for label in prediction]

        tp, fp, tn, fn = perf_measure(as_true_binary, y_pred_binary)
        aa_list.append((amino_acid, tn, fp, fn, tp))

        # Add TN, FP, FN, TP counts for each amino acid
        metrics_dict[amino_acid] = [tn, fp, fn, tp]

    # Calculate micro-average precision, recall, and F1-score
    all_tn = 0
    all_fp = 0
    all_fn = 0
    all_tp = 0
    for value_list in metrics_dict.values():
        all_tn += value_list[0]
        all_fp += value_list[1]
        all_fn += value_list[2]
        all_tp += value_list[3]

    P_microaverage_value = all_tp / (all_tp + all_fp)
    R_microaverage_value = all_tp / (all_tp + all_fn)
    F1_microaverage_value = 0
    if P_microaverage_value != 0 and R_microaverage_value != 0:
        F1_microaverage_value = 2 * (
                    (P_microaverage_value * R_microaverage_value) / (P_microaverage_value + R_microaverage_value))
    else:
        pass
    return {
        'Precision (micro)': round(P_microaverage_value, 5),
        'Recall (micro)': round(R_microaverage_value, 5),
        'F1-score (micro)': round(F1_microaverage_value, 5),
        'Per class metrics': per_class}


# Load a benchmarking dataset of A-domains and their original amino acids (chunk of ADD database)
ground_truth_data = pd.read_csv('ADD_benchmarking_dataset.tsv', sep='\t')

# Load the table containing all the predictions from SANDPUMA
sandpuma_predictions = pd.read_csv('sandpuma_predictions.tsv', sep='\t')

# Reformat original table, keeping only the selected columns that will be used for merging
# and modifying the substrate's name
selected_columns_puma = ['genomic_id', 'locus_tag', 'module', 'substrate']
origin_puma = ground_truth_data[selected_columns_puma].copy()
origin_puma['substrate'] = origin_puma['substrate'].apply(lambda x: replace_characters(x))
origin_puma.rename(columns={'substrate': 'true'}, inplace=True)
#origin_puma['genomic_id'] = origin_puma['genomic_id'].str.split('.').str[0]

# Reformat predictions table by modifying the substrate's name
sandpuma_predictions['ASM'] = sandpuma_predictions['ASM'].apply(lambda x: ast.literal_eval(x))
sandpuma_predictions['SVM'] = sandpuma_predictions['SVM'].apply(lambda x: ast.literal_eval(x))
sandpuma_predictions['pHMM'] = sandpuma_predictions['pHMM'].apply(lambda x: ast.literal_eval(x))
sandpuma_predictions['SANDPUMA'] = sandpuma_predictions['SANDPUMA'].apply(lambda x: ast.literal_eval(x))
prediction_puma = sandpuma_predictions.copy()
prediction_puma['ASM'] = prediction_puma['ASM'].apply(lambda x: [replace_characters(s) for s in x])
prediction_puma['SVM'] = prediction_puma['SVM'].apply(lambda x: [replace_characters(s) for s in x])
prediction_puma['pHMM'] = prediction_puma['pHMM'].apply(lambda x: [replace_characters(s) for s in x])
prediction_puma['SANDPUMA'] = prediction_puma['SANDPUMA'].apply(lambda x: [replace_characters(s) for s in x])

# Merge original and predicted amino acids into 1 table
df2 = origin_puma.merge(prediction_puma, on=['genomic_id', 'locus_tag', 'module'])
puma_true = df2['true'].tolist()
puma_asm = df2['ASM'].tolist()
puma_svm = df2['SVM'].tolist()
puma_hmm = df2['pHMM'].tolist()
puma_puma = df2['SANDPUMA'].tolist()

# Calculate micro-averaged evaluation metrics for each of the algorithms
metrics_asm = calculate_metrics(puma_true, puma_asm)
metrics_svm = calculate_metrics(puma_true, puma_svm)
metrics_hmm = calculate_metrics(puma_true, puma_hmm)
metrics_puma = calculate_metrics(puma_true, puma_puma)


def choose_metrics(data, metric_type):
    metric_labels = {'micro': '(micro)'}
    selected_keys = [f'{metric} {metric_labels[metric_type]}' for metric in ['Precision', 'Recall', 'F1-score']]
    prfmm_metrics = {key: value for key, value in data.items() if key in selected_keys}
    return prfmm_metrics


# Print calculated metrics into a table
stats_PRF_micro = pd.DataFrame.from_dict(
    {'SANDPUMA (Ensemble)': choose_metrics(metrics_puma, 'micro'),
     'NRPSsp (pHMM)': choose_metrics(metrics_hmm, 'micro'),
     'NRPSPredictor2 (ASM)': choose_metrics(metrics_asm, 'micro'),
     'NRPSPredictor2 (SVM)': choose_metrics(metrics_svm, 'micro')}, orient='index')

stats_PRF_micro.to_csv('prf_metrics_micro_all.tsv', sep='\t')
