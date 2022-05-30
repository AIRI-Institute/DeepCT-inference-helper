#!/usr/bin/env python

import argparse

import joblib as jl
import numpy as np


# empirical value, in range from 0 to 1; 
# the higher, the more confident prediction will be
# I recommend starting from 0.8 and increase threshold 
# if too many positive results are reported
# -- V. Fishman
THRESHOLD=0.8

# again, very empirical list
# of cell types where I saw too many significant differences
# and thus I recommend excluding them
# -- V. Fishman
EXCLUDED_CELL_TYPES = np.array([
      7,  11,  22,  23,  46,  54,  99, 117, 118, 141, 186, 192, 193,
    210, 216, 222, 226, 253, 298, 323, 334, 399, 417, 433, 445, 448,
    455, 457, 475, 500, 505, 533, 554, 574, 577, 579, 596, 613, 630,
    640, 664, 668, 669, 678, 686, 711, 716, 725, 789
])


def get_target_features_and_cell_types(distinct_features_path, target_features_path):
    def _parse_distinct_feature(distinct_feature):
            """
            Parse a combination of `cell_type|feature_name|info` into
            `(feature_name, cell_type)`
            """
            feature_description = distinct_feature.split("|")
            feature_name = feature_description[1]
            cell_type = feature_description[0]
            addon = feature_description[2]
            if addon != "None":
                cell_type = cell_type + "_" + addon
            return feature_name, cell_type


    with open(distinct_features_path) as f:
            distinct_features = list(map(lambda x: x.rstrip(), f.readlines()))

    with open(target_features_path) as f:
            target_features = list(map(lambda x: x.rstrip(), f.readlines()))

    distinct_features = [
                    i
                    for i in distinct_features
                    if _parse_distinct_feature(i)[0] in target_features
                ]

    _cell_types = []
    for distinct_feature_index, distinct_feature in enumerate(
                    distinct_features
                ):
                    feature_name, cell_type = _parse_distinct_feature(distinct_feature)
                    if feature_name not in target_features:
                        continue
                    if cell_type not in _cell_types:
                        _cell_types.append(cell_type)
    return distinct_features, target_features, _cell_types


def get_predictions(predicted_data_array, fit_logistic_regression, cell_index):
    # these are proxies for a working regulatory site in a given cell type
    ref_prediction = fit_logistic_regression.predict(predicted_data_array[:,0,cell_index,:])
    alt_prediction = fit_logistic_regression.predict(predicted_data_array[:,1,cell_index,:])
    diff_preds = (ref_prediction - alt_prediction)
    
    ref_prediction_probs = fit_logistic_regression.predict_proba(predicted_data_array[:,0,cell_index,:]).max(axis=1)
    alt_prediction_probs = fit_logistic_regression.predict_proba(predicted_data_array[:,1,cell_index,:]).max(axis=1)
    min_probs = np.minimum(ref_prediction_probs, alt_prediction_probs)

    return (diff_preds, min_probs)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Process DeepCT predictions')

    parser.add_argument('--predictions', '-p', required=True,
                        help='Path to a DeepCT concatenated_flat_predictions.npy')
    parser.add_argument('--input-vcf', '-v', required=True,
                        help='Path to an input VCF file WITHOUT THE HEADER for these predictions.')
    parser.add_argument('--output-vcf', '-o', required=True,
                        help='Path to an output VCF file. It also will be without the header.')
    
    parser.add_argument('--regression-model', default='models/logreg-score.2022-05-04.joblib', required=False,
                        help='Path to a logreg model to apply to DeepCT predictions.')
    parser.add_argument('--feature-names', default='descriptions/target_features.txt', required=False,
                        help='Path to a text file describing feature names.')
    parser.add_argument('--track-names', default='descriptions/distinct_features_nonTreated.qcfiltered.txt', required=False,
                        help='Path to a text file describing tracks (cell type * feature).')

    parser.add_argument('--include-wo-predictions', action='store_true', required=False,
                        help='Include sites without predictions in the output.')
    parser.add_argument('--report-threshold', '-t', default=THRESHOLD, type=float, required=False,
                        help='Only process predictions if their probability is higher than that.')

    return parser.parse_args()


def main():
    args = parse_arguments()

    predicted_data_array = np.load(args.predictions)
    
    fit_logistic_regression = jl.load(args.regression_model)

    distinct_features, target_features, _cell_types = get_target_features_and_cell_types(
        distinct_features_path=args.track_names,
        target_features_path=args.feature_names
    )

    var_results = [ {} for i in range(predicted_data_array.shape[0]) ]
    for cell_index in range(predicted_data_array.shape[2]):
    
        if cell_index in EXCLUDED_CELL_TYPES:
            continue
        
        (diff_preds, min_probs) = get_predictions(
            predicted_data_array, 
            fit_logistic_regression,
            cell_index
        )

        for variant in range(predicted_data_array.shape[0]):
            prob = min_probs[variant]
            if prob < args.report_threshold:
                continue

            consequence = ''
            diff = diff_preds[variant]
            if diff == 1:
                consequence = f"L{prob:.2f}" # "loss" of regulatory site
            elif diff == -1:
                consequence = f"G{prob:.2f}" # "gain"
            else:
                continue
            
            cell_type = _cell_types[cell_index]
            var_results[variant][cell_type] = consequence

    vcfout = open(args.output_vcf, 'w')
    with open(args.input_vcf,  'r') as vcfin:
        for variant in var_results:
            vcf_line = vcfin.readline()

            if variant:
                fields = vcf_line.split("\t")
                fields[7] += ";DEEPCT="+",".join( f"{variant[c]}@{c}" for c in variant.keys() )

                merged_line = "\t".join(fields)
                vcfout.write(merged_line)
            else:
                if args.include_wo_predictions:
                    vcfout.write(vcf_line)
    vcfout.close()


if __name__ == '__main__':
    main()
