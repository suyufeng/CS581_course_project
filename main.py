from Bio import SeqIO
import numpy as np
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet
import matplotlib.pyplot as plt
import os
import multiprocessing
import pandas as pd
import seaborn as sns
import scipy
import time
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo import PhyloXML
from Bio.Phylo import *
from Bio import Phylo
import sys

def remove_gap(records):
    length = len(records[0].seq)
    remain_ids = []
    for i in range(length):
        flag = False
        for record in records:
            c = record.seq[i]
            if c != "-":
                flag = True
            if flag:
                break
        if flag:
            remain_ids.append(i)
    for record in records:
        str = ""
        for id in remain_ids:
            str += record.seq[id]
        record.seq = Seq(str,SingleLetterAlphabet)
    return


def random_split(input_dir, input_file_name,
                 split_1_save_dir, split_1_file_name,
                 split_2_save_dir, split_2_file_name,
                 split_all_save_dir, split_all_file_name, tree_name = None):
    seed = 19971024
    for name in ["1000M1", "1000M4", "1000S1", "1000S4"]:
        for i in range(20):
            file_name = input_dir + "/" + name + "/R" + str(i) + "/" + input_file_name
            records = list(SeqIO.parse(file_name, "fasta"))
            fasta_names = [i.id for i in records]
            if tree_name is None:
                ids = np.arange(1000)
                np.random.seed(seed=seed)
                np.random.shuffle(ids)
            else:
                trees = list(Phylo.parse(input_dir + "/" + name + "/R" + str(i) + "/" + tree_name, 'newick'))
                term_names = [term.name for term in trees[0].get_terminals()]
                ids = [fasta_names.index(name) for name in term_names]
            first_ids = list(ids[:500])
            second_ids = list(ids[500:])
            first_records = [records[i] for i in first_ids]
            second_records = [records[i] for i in second_ids]
            remove_gap(first_records)
            remove_gap(second_records)
            SeqIO.write(first_records,split_1_save_dir + "/" + name + "/R" + str(i) + "/" + split_1_file_name, "fasta")
            SeqIO.write(second_records,split_2_save_dir + "/" + name + "/R" + str(i) + "/" + split_2_file_name, "fasta")
            first_records.extend(second_records)
            SeqIO.write(first_records,split_all_save_dir + "/" + name + "/R" + str(i) + "/" + split_all_file_name, "fasta")
            print("{} R{}".format(name, i))


def pure_sequence():
    seed = 19971024
    for name in ["1000M1", "1000M4", "1000S1", "1000S4"]:
        for i in range(20):
            file_name = "../data/" + name + "/" + "R" + str(i) + "/rose.aln.true.fasta"
            records = list(SeqIO.parse(file_name, "fasta"))
            for j in range(len(records)):
                remove_gap([records[j]])
            SeqIO.write(records, "../data/" + name + "/R" + str(i) + "/pure.fasta", "fasta")
            print("{} R{}".format(name, i))
            # exit(0)

def run_pasta():
    for name in ["1000M4"]:
        for i in [16]:
            # if  os.path.isdir("../output/pasta/" + name + "_R" + str(i)):
            #     return
            os.system("mkdir ../output/pasta/" + name + "_R" + str(i))
            os.system("python software/pasta/run_pasta.py -i ../data/" + name + "/R"+ str(i) + "/pure.fasta" +
                      " -o ../output/pasta/" + name + "_R" + str(i) + " -j test --num-cpus 10&")
            print(name + "_R" + str(i))

def run_mafft_commend(parameters):
    start_time = time.time()
    (input_dir, output_dir, tree_dir) = parameters
    if tree_dir is None:
        os.system("mafft --merge id.txt " + input_dir + " > " + output_dir)
    else:
        print("mafft --treein " + tree_dir + " --merge id_full.txt " + input_dir + " > " + output_dir )
    name = output_dir.split("/")[-2:]
    running_time = time.time() - start_time
    return (name[0] + "/" + name[1], running_time)

def run_mafft(input_dir, input_name, output_dir, tree_name = None):
    os.system("mkdir " + output_dir)
    parameters = []
    for name in ["1000M1", "1000M4", "1000S1", "1000S4"]:
        os.system("mkdir " + output_dir + "/" + name)
        for i in range(20):
            input_add = input_dir + "/" + name + "/" + "R" + str(i) + "/" + input_name
            output_add = output_dir + "/" + name + "/R" + str(i)
            tree_add = None
            if tree_name is not None:
                tree_add = input_dir + "/" + name + "/R" + str(i) + "/" + tree_name
            parameters.append((input_add, output_add, tree_add))
    pool = multiprocessing.Pool(processes=20)
    names = []
    times = []
    for name, time in pool.imap_unordered(run_mafft_commend, parameters):
        names.append(name)
        times.append(time)
    data = {}
    data["names"] = names
    data["times"] = times
    data = pd.DataFrame(data)
    data.to_csv(output_dir + "/running_time.csv", index=False)

class TwoInputModel:
    def __init__(self, input_dir, input_name1, input_name2, output_dir, num_process = 20):
        self.input_dir = input_dir
        self.input_name1 = input_name1
        self.input_name2 = input_name2
        self.output_dir = output_dir
        self.num_process = num_process

    def run_commend(self, parameters):
        start_time = time.time()
        (input_dir1, input_dir2, output_dir) = parameters
        os.system("opal --in " + input_dir1 + " --in2 " + input_dir2 + " --out " + output_dir)
        name = output_dir.split("/")[-2:]
        running_time = time.time() - start_time
        return (name[0] + "/" + name[1], running_time)

    def run(self):
        os.system("mkdir " + self.output_dir)
        parameters = []
        for name in ["1000M1", "1000M4", "1000S1", "1000S4"]:
            os.system("mkdir " + self.output_dir + "/" + name)
            for i in range(20):
                input_add1 = self.input_dir + "/" + name + "/" + "R" + str(i) + "/" + self.input_name1
                input_add2 = self.input_dir + "/" + name + "/" + "R" + str(i) + "/" + self.input_name2
                output_add = self.output_dir + "/" + name + "/R" + str(i)
                parameters.append((input_add1, input_add2, output_add))
        pool = multiprocessing.Pool(processes=self.num_process)
        names = []
        times = []
        for name, time in pool.imap_unordered(self.run_commend, parameters):
            names.append(name)
            times.append(time)
        data = {}
        data["names"] = names
        data["times"] = times
        data = pd.DataFrame(data)
        data.to_csv(self.output_dir + "/running_time.csv", index=False)

class Opal(TwoInputModel):
    def run_commend(self, parameters):
        start_time = time.time()
        (input_dir1, input_dir2, output_dir) = parameters
        os.system("opal --in " + input_dir1 + " --in2 " + input_dir2 + " --out " + output_dir)
        name = output_dir.split("/")[-2:]
        running_time = time.time() - start_time
        return (name[0] + "/" + name[1], running_time)

class Muscle(TwoInputModel):
    def run_commend(self, parameters):
        start_time = time.time()
        (input_dir1, input_dir2, output_dir) = parameters
        os.system("muscle -profile -in1 " + input_dir1 + " -in2 " + input_dir2 + " -out " + output_dir)
        name = output_dir.split("/")[-2:]
        running_time = time.time() - start_time
        return (name[0] + "/" + name[1], running_time)

class Clustalw(TwoInputModel):
    def run_commend(self, parameters):
        start_time = time.time()
        (input_dir1, input_dir2, output_dir) = parameters
        os.system("clustalw2 -PROFILE -PROFILE1=" + input_dir1 + " -PROFILE2=" + input_dir2 + " -OUTFILE=" + output_dir + " -OUTPUT=fasta")
        name = output_dir.split("/")[-2:]
        running_time = time.time() - start_time
        return (name[0] + "/" + name[1], running_time)

class GuideTree:
    def __init__(self, input_dir, input_name, output_dir, tree_dir, num_process = 20):
        self.input_dir = input_dir
        self.input_name = input_name
        self.output_dir = output_dir
        self.tree_dir = tree_dir
        self.num_process = num_process

    def run_commend(self, parameters):
        start_time = time.time()
        (input_dir1, input_dir2, output_dir) = parameters
        os.system("opal --in " + input_dir1 + " --in2 " + input_dir2 + " --out " + output_dir)
        name = output_dir.split("/")[-2:]
        running_time = time.time() - start_time
        return (name[0] + "/" + name[1], running_time)

    def run(self):
        parameters = []
        for name in ["1000M1", "1000M4", "1000S1", "1000S4"]:
            os.system("mkdir " + self.output_dir + "/" + name)
            for i in range(20):
                input_add = self.input_dir + "/" + name + "/" + "R" + str(i) + "/" + self.input_name
                tree = self.tree_dir + " "
                output_add = self.output_dir + "/" + name + "/R" + str(i)
                parameters.append((input_add, output_add))
        pool = multiprocessing.Pool(processes=self.num_process)
        names = []
        times = []
        for name, time in pool.imap_unordered(self.run_commend, parameters):
            names.append(name)
            times.append(time)
        data = {}
        data["names"] = names
        data["times"] = times
        data = pd.DataFrame(data)
        data.to_csv(self.output_dir + "/running_time.csv", index=False)

def mv_files():
    for name in ["1000M1", "1000M4", "1000S1", "1000S4"]:
        for i in range(20):
            os.system("cp ../output/pasta/" + name + "_R" + str(i) + "/test.tre ../data/" + name + "/R" + str(i) + "/pasta_tree")

def calculate(recordsa, recordsb):
    occur_ids = []
    keys = []
    for key in recordsa.keys():
        ids = []
        keys.append(key)
        for index, seq in enumerate(recordsa[key].seq):
            if seq != '-':
                ids.append(index)
        occur_ids.append(ids)

    strs = []
    for key in keys:
        seq_key = ""
        for str in recordsb[key].seq:
            seq_key += str
        strs.append(seq_key)
    num = [-1] * len(strs)
    cnt = 0
    correct_cnt = 0
    for row_id in range(len(strs[0])):
        charactors = []
        for index, str in enumerate(strs):
            if str[row_id] != '-':
                num[index] += 1
                charactors.append((index, num[index]))
        for index1, num1 in charactors:
            for index2, num2 in charactors:
                if index1 == index2:
                    continue
                cnt += 1
                if occur_ids[index1][num1] == occur_ids[index2][num2]:
                    correct_cnt += 1
    return 1. * correct_cnt / cnt

def get_metrics(true_dir, aligned_dir, id, name):
    true_records = SeqIO.to_dict(SeqIO.parse(true_dir, "fasta"))
    aligned_records = SeqIO.to_dict(SeqIO.parse(aligned_dir, "fasta"))

    key_test = 0
    for _ in true_records.keys():
        key_test = _
        break

    compression_rate = 1. * len(aligned_records[key_test].seq) / len(true_records[key_test].seq)

    times = []
    times.append(time.time())
    spfp = calculate(aligned_records, true_records)
    times.append(time.time())
    distance = get_kl_dis(true_dir, aligned_dir)
    times.append(time.time())
    spfn = calculate(true_records, aligned_records)


    return (spfp, spfn, compression_rate, distance, times[2] - times[1], times[1] - times[0], id, name)

def get_gaps(records, keys):
    results = []
    for key in keys:
        record = records[key]
        poss = []
        for index, seq in enumerate(record.seq):
            if seq != '-':
                poss.append(index)
        for index, pos in enumerate(poss):
            if index == 0:
                results.append(pos + 1)
            else:
                results.append(pos - poss[index - 1])
        results.append(len(record) - poss[-1])

    return np.array(results)

def KL(P,Q):
     epsilon = 0.00001
     P = P-1+epsilon
     Q = Q-1+epsilon
     P = P /sum(P)
     Q = Q / sum(Q)
     divergence = np.sum(P*np.log(P/Q))
     print(divergence)
     return divergence

def get_kl_dis(true_dir, aligned_dir, id=None, name=None):
    true_records = SeqIO.to_dict(SeqIO.parse(true_dir, "fasta"))
    aligned_records = SeqIO.to_dict(SeqIO.parse(aligned_dir, "fasta"))
    keys = [key for key in true_records.keys()]
    true_gaps = get_gaps(true_records, keys)
    aligned_gaps = get_gaps(aligned_records, keys)
    if id is None:
        return KL(aligned_gaps, true_gaps)
    else:
        return KL(aligned_gaps, true_gaps), id, name

def get_summary_commend(parameters):
    (name, method, type, i) = parameters
    return get_metrics("../data/" + name + "/R" + str(i) + "/rose.aln.true.fasta", "../output/"
                                                                     + method + "_" + type + "/" + name + "/R" + str(i), i, name)
def get_summary():
    for method in ["clustal", "mafft", "muscle", "opal"]:
        for type in ["est_tree", "random_pasta", "random_true", "true_tree"]:
    # for method in ["muscle"]:
    #     for type in ["random_true"]:
            print("finished {}/{}".format(method, type))
            data = {}
            table = pd.read_csv("../output/" + method + "_" + type + "/running_time.csv")
            names = table['names']
            times = table['times']
            name2time = {}
            for index, name in enumerate(names):
                name2time[name] = times[index]

            names = []
            times = []
            spfps = []
            spfns = []
            compressions = []
            distances = []
            time_distances = []
            time_spfns = []

            parameters = []
            for name in ["1000M1", "1000M4", "1000S1", "1000S4"]:
                for i in range(20):
                    parameters.append((name, method, type, i))
            # get_summary_commend(parameters[0])

            pool = multiprocessing.Pool(processes=20)

            for spfp, spfn, compression, distance, time_spfn, time_distance, i, name in pool.imap_unordered(get_summary_commend, parameters):
                names.append(name + "/R" + str(i))
                times.append(name2time[name + "/R" + str(i)])
                spfps.append(spfp)
                spfns.append(spfn)
                distances.append(distance)
                time_distances.append(time_distance)
                time_spfns.append(time_spfn)
                compressions.append(compression)

            data["name"] = names
            data["running_time"] = np.array(times).round(decimals=5)
            data["spfp"] = np.array(spfps).round(decimals=5)
            data["spfn"] = np.array(spfns).round(decimals=5)
            data["compression"] = np.array(compressions).round(decimals=5)
            data["distance"] = np.array(distances).round(decimals=5)
            data["time_distance"] = np.array(time_distances).round(decimals=5)
            data["time_spfn"] = np.array(time_spfns).round(decimals=5)
            data = pd.DataFrame(data)
            data.sort_values(by='name').reset_index()
            data.to_csv("../result/" + method + "_" + type + ".csv", index=False, columns=["name", "running_time", "spfp", "spfn", "compression", "distance", "time_distance", "time_spfn"])

def get_kl_dis_commend(parameters):
    (name, method, type, i) = parameters
    return get_kl_dis("../data/" + name + "/R" + str(i) + "/rose.aln.true.fasta", "../output/"
                                                                     + method + "_" + type + "/" + name + "/R" + str(i), i, name)

def change_kl_dis():
    for method in ["clustal", "mafft", "muscle", "opal"]:
        for type in ["est_tree", "random_pasta", "random_true", "true_tree"]:
            data = pd.read_csv("../result/" + method + "_" + type + ".csv")
            kl_divs = [0] * 80
            names = list(data['name'].to_numpy())
            parameters = []
            for name in ["1000M1", "1000M4", "1000S1", "1000S4"]:
                for i in range(20):
                    parameters.append((name, method, type, i))

            pool = multiprocessing.Pool(processes=20)

            for kl_div, i, name in pool.imap_unordered(get_kl_dis_commend, parameters):
                kl_divs[names.index(name + "/R" + str(i))] = kl_div

            tmp = data['time_distance'].to_numpy().copy()
            data['time_distance'] = data['time_spfn'].to_numpy().copy()
            data['time_spfn'] = tmp
            data['distance'] = kl_divs
            data['compression'] = 1. / data['compression'].to_numpy()
            for col in data.columns:
                if col == 'name':
                    continue
                data[col] = np.round(data[col].to_numpy(), decimals=4)
            data.to_csv("../new_result/" + method + "_" + type + ".csv", index=False,
                        columns=["name", "running_time", "spfp", "spfn", "compression", "distance", "time_distance",
                                 "time_spfn"], sep='\t')

def get_key_value(names, name, key, data, mean_flag = True):
    result = []
    for i in range(20):
        id = names.index(name + "/R" + str(i))
        result.append(data[key][id])
    if mean_flag:
        return np.round(np.average(np.array(result)), decimals=4)
    else:
        return np.round(np.array(result), decimals=4)

def print_table(keys):
    clustal_spfn = []
    mafft_spfn = []
    # for method in ["clustal", "mafft", "muscle", "opal"]:
    t1 = []
    t2 = []
    for method in [ "mafft", "muscle"]:
    # for method in ["muscle"]:
        strr = method + "&"
        ts = []
        # for type in ["true_tree", "random_true", "random_pasta", "est_tree"]:

        for type in ["est_tree"]:
            data = pd.read_csv("../new_result/" + method + "_" + type + ".csv", sep ='\t')
            names = list(data['name'].to_numpy())
            # for name in ["1000S4", "1000M4", "1000S1", "1000M1"]:

            for name in ["1000M1"]:
                for key in keys:
                    t = get_key_value(names, name, key, data, mean_flag=False)
                    ts.extend(t)
                    strr += "{}&".format(t)
        if method == 'mafft':
            t1 = ts.copy()
        else:
            t2 = ts.copy()
        ts = []

        strr = strr[:-1] + "\\\\"
        print(strr)
    draw(t1, t2, "../figure/muscle_mafft_distance")

def draw(our, RCK, output_dir):
    # %%

    # independent two samples t-test
    xx = scipy.stats.ttest_ind(our, RCK)
    print(scipy.stats.ttest_ind(our, RCK))
    print(xx[1] / 2)

    # wilcoxon signed-rank test
    xx = scipy.stats.wilcoxon(our, RCK)
    print(scipy.stats.wilcoxon(our, RCK))
    print(xx[1] / 2)

    # sns.set(style="whitegrid")
    sns.set_style("whitegrid", {'axes.grid': False})
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    fig, ax = plt.subplots(figsize=(10, 10))
    grid = plt.GridSpec(1, 1)

    ax = plt.subplot(grid[0, 0])
    sns.scatterplot(x=RCK, y=our, ax=ax, s=120, alpha=0.7)
    # ax.set_title("Pearson correlation of predicted and measured intensities", fontdict={"fontsize": 40}, pad=30)
    ax.set_xlabel("Mafft-merge", fontdict={"fontsize": 20})
    ax.set_ylabel("muscle", fontdict={"fontsize": 20})
    # ax.set_xticklabels(fontdict={"fontsize": 20})

    ax.set_xlim(7.80, 10.1)
    ax.set_ylim(7.80, 10.1)

    ax.plot((0, 20), (0, 20), "--", color="lightgrey")

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(15)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(15)

    props = dict(boxstyle='round', facecolor='gray', alpha=0.15)
    textstr = "Mafft-merge > muscle:\n" + r"$p = 0.3544$"
    ax.text(0.58, 0.15, textstr, transform=ax.transAxes, fontsize=17,
            verticalalignment='top', bbox=props)

    plt.savefig(output_dir + ".pdf")

# ["name", "running_time", "spfp", "spfn", "compression", "distance", "time_distance",
#                                  "time_spfn"]
print_table(['distance'])

# change_kl_dis()

# random_split("../data", "pasta_estimated.fasta", "../data", "pasta_split_0.fasta", "../data", "pasta_split_2.fasta", "../data", "pasta_split_all.fasta")

# random_split("../data", "pasta_estimated.fasta", "../data", "pasta_split_est_tree_0.fasta", "../data", "pasta_split_est_tree_1.fasta", "../data", "pasta_split_est_tree_all.fasta", tree_name="pasta_tree")
# random_split("../data", "rose.aln.true.fasta", "../data", "true_tree_split_0.fasta", "../data", "true_tree_split_1.fasta", "../data", "true_tree_split_all.fasta", tree_name="rose.mt")

# change_tree_form("../data", "pasta_estimated.fasta","pasta_tree", "")
# run_mafft("../data", "pasta_split_est_tree_all.fasta", "../output/mafft_est_tree")
# run_pasta()
# run_opal("../data", "random_split_0.fasta", "random_split_1.fasta", "../output/opal_random_true")
# a = Muscle("../data", "pasta_split_est_tree_0.fasta", "pasta_split_est_tree_1.fasta", "../output/muscle_est_tree")

# a = Muscle("../data", "true_tree_split_0.fasta", "true_tree_split_1.fasta", "../output/muscle_true_tree")
# a = Clustalw("../data", "pasta_split_1.fasta", "pasta_split_2.fasta", "../output/clustal_random_pasta")

# a = Muscle("../data", "random_split_0.fasta", "random_split_1.fasta", "../output/muscle_random_true")
# a.run()

# get_metrics("../data/1000M1/R0/rose.aln.true.fasta", "../output/clustal_est_tree/1000M1/R0")

# change_kl_dis()