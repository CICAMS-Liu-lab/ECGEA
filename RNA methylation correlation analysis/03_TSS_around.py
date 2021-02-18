import matplotlib
matplotlib.use('Agg')
import sys
sys.path.append("/usr/local/python3/lib/python3.5/site-packages")
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import glob
import argparse
import matplotlib as mpl
from my_functions import slide_window

"""
基于滑窗法
从CX_report.txt.simplified中提取出来特定基因的信息
"""

def parse_args():
    parser = argparse.ArgumentParser(description="methylation level around TSS")
    parser.add_argument("-wd",default="/root/mnt/analysis/baidu_II/WGBS_data",help="CX_report.txt.simplified")
    parser.add_argument("-sd",default="/root/mnt/analysis/baidu_II/RNAseq/S100A14",help="save dir")
    parser.add_argument("-tl",help="tumor list")
    parser.add_argument("-nl", help="normal list")
    parser.add_argument("-sn", help="save name")
    parser.add_argument("-w", type=int, default=500,help="window length")
    parser.add_argument("-b", type=int, default=130, help="bins")
    args = parser.parse_args()
    return args.wd,args.sd,args.tl,args.nl,args.sn,args.w,args.b

def process(samples,wd,wl,bins):
    methylation_level = []
    for sample in samples:
        try:
            df = pd.read_table("%s/%s/%s.CX_report.txt.simplified" %(wd,sample,sample),header=None,names=["chr","loc","meth","unmeth","pvalue"],index_col=0)
            result = slide_window("chr1", 153586731, 153589462, "-", df, bin_length=wl,bins=bins)
            methylation_level.append(result)
        except:
            pass
    df = pd.DataFrame(methylation_level)
    mean_value = df.mean()
    return mean_value

def plot(mean1,mean2,sd,sn,wl,bins):
    mpl.rcParams["font.size"] = 16
    mpl.rcParams["font.weight"] = "bold"
    f, ax = plt.subplots(figsize=(10, 10))
    ax.plot(range(0, bins), mean1, "red", label="Tumor")
    ax.plot(range(0, bins), mean2, "blue", label="Normal")
    raw_xticks = np.array([0,15,35,55,75,95,115,130])
    modify_xticks = np.ceil(raw_xticks*bins/130)
    ax.set_xticks(modify_xticks.tolist())
    ax.set_xticklabels([-15,0,20,40,60,80,100,115])
    ax.set_xlabel("Distance from TSS(% gene length)", fontdict={"size": 20, "weight": "bold"})
    ax.set_ylabel("Methylation level", fontdict={"size": 20, "weight": "bold"})
    ax.legend()
    f.savefig("%s/%s_%s_%s_tss.png" %(sd,sn,wl,bins),dpi=100,bbox_inches='tight')

def pipeline(wd,sd,tl,nl,sn,wl,bins):
#    os.chdir(wd)
    tumorList = [i.rstrip() for i in open(tl)]
    normalList = [i.rstrip() for i in open(nl)]
#    tumor_mean_value = process(tumorList,wl,bins)
#    normal_mean_value = process(normalList,wl,bins)
    tumor_mean_value = process(tumorList,wd,wl,bins)
    normal_mean_value = process(normalList,wd,wl,bins)
    plot(tumor_mean_value, normal_mean_value, sd, sn, wl, bins)

if __name__ == "__main__":
    wd,sd,tl,nl,sn,wl,bins = parse_args()
    pipeline(wd,sd,tl,nl,sn,wl,bins)
