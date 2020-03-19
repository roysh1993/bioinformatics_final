from settings import *
import json
import os
import shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def generate_plots(stats_list, directory = OUTPUT_DIR):
    cov = sorted(stats_list.keys())
    N50_scaffolds =list( map(lambda x:np.mean([n["Scaffold Stats"]["N50"] for n in stats_list[x]]), cov))
    N50_contigs = list(map(lambda x:np.mean([n["Contig Stats"]["N50"] for n in stats_list[x]]), cov))

    print(N50_scaffolds)
    # scaffolds
    plt.plot(cov, N50_scaffolds)
    plt.xlabel("coverage")
    plt.ylabel("N50 avg")
    plt.title("N50 per cov")
    plt.savefig(os.path.join(directory,"N50_scaff.png"))
    #contigs
    plt.plot(cov, N50_contigs)
    plt.xlabel("coverage")
    plt.ylabel("N50 avg")
    plt.title("N50 per cov")
    plt.savefig(os.path.join(directory,"N50_contig.png"))
    return

def save_stats(stats,directory_name):
    file_path= os.path.join(directory_name,"stats.json")
    if not os.path.exists(directory_name):
        os.makedirs(directory_name)
    with open(file_path, 'w') as fp:
        json.dump(stats, fp)
    return


# def temp_func():
#     rootdir = '/home/roy.shafir/res/Staphylococcus_aureus_output'
#     # rootdir = "/Users/royshafir/Google Drive/CS _Degree/Technion/First Year/Intro To Bioinformatics/Final Project/bioinformatics_final"
#     for subdir, dirs, files in os.walk(rootdir):
#         if subdir.endswith(".git"):
#             print("HEREEEE")
#             print(subdir)
#             continue
#
#         if subdir.endswith("corrected"):
#             shutil.rmtree(subdir)
#             print("HERE")
#         # for file in files:
#         #     a = os.path.join(subdir, file)
#         #     if not file.endswith(".fasta") or file.endswith(".json"):
#         #         os.remove(a)
#
# if __name__ =='__main__':
#     temp_func()
