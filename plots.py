from settings import *
import json
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def generate_plots(stats_list,N50_list, directory = OUTPUT_DIR):
    cov = sorted(stats_list.keys())
    N50_scaffolds = map(lambda x:np.mean([n["Scaffold Stats"]["N50"] for n in stats_list[x]]), cov)
    N50_contigs = map(lambda x:np.mean([n["Contig Stats"]["N50"] for n in stats_list[x]]), cov)

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
