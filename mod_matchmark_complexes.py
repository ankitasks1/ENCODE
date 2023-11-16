# Read in known pairs of stringDB data in encode data and calculate precison and recall (filtered for the dataset)

import sys

str = sys.argv[1]
ens = sys.argv[2]
values = [0, 50, 100, 150, 200, 250, 300, 323]

def predictor(ens, str, values):
    for val in values:
        # print(val)
        # Read in encode data 
        ens_dict = dict()
        

        with open(ens) as e:
            for x in e.readlines():
                x = x.rstrip()
                v = x.split('\t')
                # print(v)
                if v[1] + '_' + v[0] in ens_dict:
                    continue
                elif float(v[2]) >= float(val):
                        ens_dict[v[0] + '_' + v[1]] = v[2]
        # print(val, len(set(ens_dict)))

        # Read in stringDB data 
        str_dict = dict()

        with open(str) as s:
            for x in s.readlines():
                x = x.rstrip()
                v = x.split('\t')
                # print(count, v)
                key1 = v[0] + '_' + v[1]
                key2 = v[1] + '_' + v[0]
                if key1 not in str_dict and key2 not in str_dict:
                    str_dict[key1] = v[2]
                        


        # if there is no elements from encode exist in stringdb (at very high score of encode) i.e. everything is novel interaction, then str_set will be empty
        # in that case TP, FP, FN cant be calculated and I will just put NA
        # print(ens_dict)
        # print(str_dict)

        ens_set = set(ens_dict)
        str_set = set(str_dict)
        if len(str_set) >= 1:
            if len(ens_set.intersection(str_set)) >= 1:
                # print('encode_set',len(ens_set))
                # print('str_set',len(str_set))
                # print('interaction_set',len(ens_set.intersection(str_set)))

                # define parameters /check
                true_positives = len(ens_set.intersection(str_set))
                false_positives = len(ens_set) - true_positives
                false_negatives = len(str_set) - true_positives
                precison = true_positives/ (true_positives + false_positives)
                recall = true_positives/ (true_positives + false_negatives)
                f1score = 2 * (precison * recall) / (precison + recall)
                # print(true_positives, false_positives, false_negatives)
                print(val, '\t',true_positives, '\t',false_positives, '\t',false_negatives,'\t',round(precison,3), '\t',round(recall,3), '\t', round(f1score, 3))
                # print('\n')
            elif len(ens_set.intersection(str_set)) < 1:
                print(val, '\t', "NA", '\t', "NA",'\t', "NA",'\t', "NA",'\t', "NA", '\t', "NA") # where NA can indicate a novel pair at that co-association score, (a score at which not even a single known pair was identified in stringDB)
        else:
            print(val, '\t', "NA",'\t', "NA",'\t', "NA",'\t', "NA", '\t', "NA", '\t', "NA") # where NA can indicate a novel pair at that co-association score, (a score at which not even a single known pair was identified in stringDB)



predictor(ens, str, values)