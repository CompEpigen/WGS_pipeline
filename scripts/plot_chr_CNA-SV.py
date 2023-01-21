import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import gridspec



parser = argparse.ArgumentParser()
parser.add_argument('--cna', type = str, help='Copy number file')
parser.add_argument('--sv', type = str, help='SV file')
parser.add_argument('--foldback', type = str, help='tsv file containing foldback inversions (format: chr pos orientation)')
parser.add_argument('--sample', type = str, help='sample name')
parser.add_argument('--sex', type = str,default="XX", help='Sex: XX or XY')
parser.add_argument('-o', type = str, help='Output basename')
parser.add_argument('--format', type = str, help='Output format: png or svg')
args = parser.parse_args()


sample = args.sample
use_foldbacks = args.foldback is not None

def add_half_ellipse(ax,x1,x2,ymax,color):
    e = patches.Ellipse(xy=((x1+x2)/2,0),width=abs(x2-x1),height=ymax*2,fill=False,edgecolor=color)
    ax.add_patch(e)

def chr_to_int(chr):
    chr=str(chr)
    if chr.isdigit():
        return int(chr)
    elif chr=="X":
        return 23
    else:
        return 24

matplotlib.rcParams.update({'font.size': 20})
matplotlib.rcParams['font.sans-serif'] = "helvetica"
matplotlib.rcParams['font.family'] = "sans-serif"
df_CNA = pd.read_csv(args.cna,sep="\t",header=None,names=["chr","start","end","ratio","color"],dtype={"chr":str})
df_SV = pd.read_csv(args.sv,sep="\t",header=0,dtype={"chr":str})
if use_foldbacks:
    df_foldback = pd.read_csv(args.foldback,sep="\t",header=None,names=["chr","pos","orientation"],dtype={"chr":str})


chromosomes = [str(x) for x in range(1,23)] + ["X"]
if "Y" in args.sex:
    chromosomes = chromosomes + ["Y"]
else:
    open(args.o+"_chrY."+args.format,"w").close() # create empty chrY file so that the pipeline knows that all output files are there.

for chr in chromosomes:
    df_CNA_chr = df_CNA.loc[df_CNA["chr"]==chr,:]
    df_SV_chr = df_SV.loc[df_SV["chr"]==chr,:]
    df_SV_chr.reset_index(inplace=True)
    if use_foldbacks:
        df_foldback_chr = df_foldback.loc[df_foldback["chr"]==chr,:]
        df_foldback_chr.reset_index(inplace=True,drop="index")

    # Set axes
    fig = plt.figure(figsize=(20,10))
    spec = gridspec.GridSpec(ncols=1, nrows=2,hspace=0.0, height_ratios=[1, 3])
    ax=[]
    ax.append(fig.add_subplot(spec[0]))
    ax.append(fig.add_subplot(spec[1]))

    ymax = np.max([2.2,np.max(df_CNA_chr["ratio"])*2*1.1])
    ymin = np.min([1.2,np.min(df_CNA_chr["ratio"])*2*0.9])
    xmax = np.max(df_CNA_chr["end"])*1.02
    print(xmax)
    

    ax[0].set_xlim(0,xmax)
    ax[0].set_ylim(0,2)
    ax[0].tick_params(axis='y', which='both', left=False,right=False,labelleft=False)
    ax[0].tick_params(axis='x', which='both', bottom=False,top=False,labelbottom=False)
    ax[0].set_ylabel("Breakpoints")
    ax[0].set_title(sample + " (chr"+str(chr)+")")


    ax[1].set_ylim(ymin,ymax)
    ax[1].set_xlim(0,xmax)
    ax[1].grid(visible=True,which="major",axis="x",color="darkgray")
    ax[1].grid(visible=True,which="minor",axis="x",color="lightgray",linestyle="dashed")
    ax[1].grid(visible=True,which="major",axis="y",color="lightgray",linestyle="dashed")
    ax[1].set_axisbelow(True)
    ax[1].set_xlabel("Position on chr"+str(chr) + " (Mb)")
    ax[1].set_ylabel("Copy number")



    chr_length = np.max(df_CNA_chr["end"])
    if chr_length>120000000:
        ax[1].xaxis.set_major_locator(plt.MultipleLocator(20000000))
        ax[1].xaxis.set_minor_locator(plt.MultipleLocator(2000000))
    else:
        ax[1].xaxis.set_major_locator(plt.MultipleLocator(10000000))
        ax[1].xaxis.set_minor_locator(plt.MultipleLocator(1000000))
    if ymax<10:
        ax[1].yaxis.set_major_locator(plt.MultipleLocator(1))
        ax[1].yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    else:
        ax[1].yaxis.set_major_locator(plt.MultipleLocator(5))
        ax[1].yaxis.set_minor_locator(plt.MultipleLocator(1))
    def format_func(value, tick_number):
            return("{:.0f}".format(value/1000000))

    def format_func2(value, tick_number):
        if value%50000000==0:
            return("{:.0f}".format(value/1000000))
        else:
            return ""

    ax[1].ticklabel_format(axis="x", style="sci", scilimits=(6,6))
    ax[1].xaxis.set_major_formatter(plt.FuncFormatter(format_func))

    # CNA plot
    ax[1].scatter(df_CNA_chr["start"],2*df_CNA_chr["ratio"],c=df_CNA_chr["color"],s=2.0,marker="o")


    # SV 
    max_SV_size=1000
    for i in range(df_SV_chr.shape[0]):
        if df_SV_chr.loc[i,"end"] - df_SV_chr.loc[i,"start"] > max_SV_size:
            max_SV_size = df_SV_chr.loc[i,"end"] - df_SV_chr.loc[i,"start"]

    for i in range(df_SV_chr.shape[0]):
        print(i)
        height = 0.2 + 1.3 * (df_SV_chr.loc[i,"end"] - df_SV_chr.loc[i,"start"]) / max_SV_size
        if df_SV_chr.loc[i,"type"]=="BND":
            top = 0.2 + 1.6 * (chr_to_int(df_SV_chr.loc[i,"chr2"]))/27
            ax[0].vlines(x=df_SV_chr.loc[i,"start"],ymin=0,ymax=top,color="green")
            ax[0].text(x=df_SV_chr.loc[i,"start"],y=top,s=df_SV_chr.loc[i,"chr2"],horizontalalignment="center",color="black",fontsize=14)
        elif df_SV_chr.loc[i,"type"]=="INV":
            add_half_ellipse(ax[0],df_SV_chr.loc[i,"start"],df_SV_chr.loc[i,"end"],height,color="black")
        elif df_SV_chr.loc[i,"type"]=="DEL":
            add_half_ellipse(ax[0],df_SV_chr.loc[i,"start"],df_SV_chr.loc[i,"end"],height,color="blue")
        elif df_SV_chr.loc[i,"type"]=="DUP":
            add_half_ellipse(ax[0],df_SV_chr.loc[i,"start"],df_SV_chr.loc[i,"end"],height,color="red")
    # Foldback inv
    if use_foldbacks:
        for i in range(df_foldback_chr.shape[0]):
            ax[0].vlines(x=df_foldback_chr.loc[i,"pos"],ymin=0,ymax=1,color="purple")
            if df_foldback_chr.loc[i,"orientation"]=="left":
                ax[0].arrow(x=df_foldback_chr.loc[i,"pos"],y=0.5,dx=-0.008*xmax,dy=0,color="purple",width=0.01,head_width=0.07,head_length=0.005*xmax,length_includes_head=True)
            else:
                ax[0].arrow(x=df_foldback_chr.loc[i,"pos"],y=0.5,dx=+0.008*xmax,dy=0,color="purple",width=0.01,head_width=0.07,head_length=0.005*xmax,length_includes_head=True)

    fig.align_ylabels(ax)
    fig.savefig(args.o+"_chr"+chr+"."+args.format,bbox_inches="tight",pad_inches=0.1,dpi=200)
    plt.cla() 
    plt.clf() 
    plt.close('all')
    #plt.show()






