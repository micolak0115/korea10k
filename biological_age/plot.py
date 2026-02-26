import os
from collections import Counter
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import PercentFormatter

def draw_barplot(group, list_age_group, outdir, figname, width=7, color="grey"):
    dict_cnt_age_group = dict(sorted(dict(Counter(list_age_group)).items()))
    list_keys = list(dict_cnt_age_group.keys())
    list_vals = list(dict_cnt_age_group.values())
    sum_vals = sum(list_vals)
    list_props = list(map(lambda x: int(x)/int(sum_vals)*100, list_vals))
    list_props_rounded = list(map(lambda x: round(x, 1), list_props))
    _, ax = plt.subplots(figsize=(5,5),layout="constrained")
    ax.bar(list_keys, list_props_rounded, width=width, color=color, zorder=111)
    for i in ax.containers:
        ax.bar_label(i, size=13)
    plt.xlabel("Age (yrs)", fontsize=18)
    plt.ylabel("Sample Proportion (%)", fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.grid(visible=True, zorder=-999)
    plt.title(group)
    figpath = os.path.join(outdir, f"age_distribution_{figname}.png")
    plt.xlim(0, 100)
    
    plt.savefig(figpath)
    plt.show()
    plt.close()

def draw_histogram_group(list_age, color="grey"):
    plt.hist(list_age, color=color, weights=np.ones(len(list_age)) / len(list_age), zorder=111)
    plt.xlabel("Age (yrs)", fontsize=18)
    plt.ylabel("Sample Proportion (%)", fontsize=18)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.grid(visible=True, zorder=-999)
    plt.show()
    plt.close()

def draw_histogram_train_test(y_train, y_test, outdir, **kwargs):
    from matplotlib.ticker import PercentFormatter
    fig, axes = plt.subplots(1,2, sharey=True)
    axes[0].hist(y_train, label="Train", color="skyblue", weights=np.ones(len(y_train)) / len(y_train), **kwargs)
    axes[1].hist(y_test, label="Test", color="orange", weights=np.ones(len(y_test)) / len(y_test), **kwargs)
    axes[0].legend(loc="best")
    axes[1].legend(loc="best")
    axes[0].set_xlim(19, 81)
    axes[1].set_xlim(19, 81)
    fig.supylabel("Proportion", fontsize=12)
    fig.supxlabel("Sample Age", fontsize=12)
    fig.suptitle("Histograms of Stratified Split Between Train and Test by Age")
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.savefig(os.path.join(outdir, "distribution_train_test_split.tiff"), dpi=600)
    plt.show()
    plt.close()