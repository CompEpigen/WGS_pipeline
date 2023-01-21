import argparse
import pycircos
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from Bio import SeqIO 
Garc    = pycircos.Garc
Gcircle = pycircos.Gcircle

parser = argparse.ArgumentParser()
parser.add_argument('--cna', type = str, help='Copy number file')
parser.add_argument('--sv', type = str, help='SV file')
parser.add_argument('--sample', type = str, help='sample name')
parser.add_argument('--sex', type = str,default="XX", help='Sex: XX or XY')
parser.add_argument('-o', type = str, help='Output file')
args = parser.parse_args()



data_dir = "data"

#Set chromosomes
circle = Gcircle() 
chr_lengths={}
with open(data_dir+"/circos_chromosome_general.csv") as f:
    f.readline()
    for line in f:
        line   = line.rstrip().split(",") 
        name   = line[0]
        length = int(line[-1]) 
        chr_lengths[name] = length
        arc    = Garc(arc_id=name, size=length, interspace=1, raxis_range=(750,800), labelposition=60, label_visible=True) # ,facecolor="white"
        circle.add_garc(arc) 

circle.set_garcs(00,352) 


#cytoband
import collections
color_dict   = {"gneg":"#FFFFFF00", "gpos25":"#EEEEEE", "gpos50":"#BBBBBB", "gpos75":"#777777", "gpos100":"#000000", "gvar":"#FFFFFF00", "stalk":"#C01E27", 
               "acen":"#D82322"}

arcdata_dict = collections.defaultdict(dict)
with open(data_dir+"/circos_chromosome_cytoband.csv") as f:
    f.readline()
    for line in f:
        line  = line.rstrip().split(",")
        name  = line[0]
        start = int(line[1])-1 
        width = int(line[2])-(int(line[1])-1) 
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["widths"]    = [] 
            arcdata_dict[name]["colors"]    = [] 
        arcdata_dict[name]["positions"].append(start) 
        arcdata_dict[name]["widths"].append(width)
        arcdata_dict[name]["colors"].append(color_dict[line[-1]])

for key in arcdata_dict:
    circle.barplot(key, data=[1]*len(arcdata_dict[key]["positions"]), positions=arcdata_dict[key]["positions"], 
                   width=arcdata_dict[key]["widths"], raxis_range=[750,800], facecolor=arcdata_dict[key]["colors"])    


color_map={"green":"#27ae60","blue":"#4a69bd","red":"#e55039","black":"black"}

#scatter plot
values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open(args.cna,"r") as f:
    f.readline()
    nskipped=0
    groupsize=0
    last_color="no"
    last_name="-1"
    for line in f:
        line  = line.rstrip().split("\t")
        name  = "chr"+line[0]     
        start = int(line[1])-1
        end   = int(line[2]) 
        mid   = (start+end)/2
        if args.sex =="XY" and name in ["chrX","chrY"]:
            value = min(float(line[3]),3.5)
        else:
            value = min(float(line[3])*2,3.5)
        values_all.append(value)
        color = line[-1].rstrip("\n") 
        if name not in arcdata_dict:
            arcdata_dict[name]["positions"] = []
            arcdata_dict[name]["values"] = []
            arcdata_dict[name]["colors"] = []
        if last_name ==last_name and last_color ==color and groupsize>10 and nskipped<=20: # skip some points, when a large group is similar...
            nskipped+=1
            groupsize+=1
        else:
            arcdata_dict[name]["positions"].append(mid) 
            arcdata_dict[name]["values"].append(value)
            arcdata_dict[name]["colors"].append(color_map[color])
            if last_name ==name and last_color==color:
                groupsize+=1
                nskipped=0
            else:
                groupsize=0
            last_name=name
            last_color=color
    
vmin, vmax = min(values_all), max(values_all) 
vmin2=vmin-0.05*abs(vmin)
vmax2= vmax+0.05*abs(vmax)

for chr in chr_lengths:
    for y in [1,2,3]:
        circle.lineplot(chr,data=[y]*20,positions=np.linspace(0,chr_lengths[chr],20),rlim=[vmin2, vmax2], raxis_range=(500,740),linecolor="grey",linewidth=0.5)

for key in arcdata_dict:
    circle.scatterplot(key, data=arcdata_dict[key]["values"], positions=arcdata_dict[key]["positions"], 
                       rlim=[vmin2,vmax2], raxis_range=(500,740), facecolor=arcdata_dict[key]["colors"], spine=True) #scatter plot


#linkplot
values_all   = [] 
arcdata_dict = collections.defaultdict(dict)
with open(args.sv,"r") as f:
    for line in f:
        line  = line.rstrip().split("\t")
        name1  = "chr"+line[0]     
        start1 = int(line[1])-1
        end1   = int(line[2])
        name2  = "chr"+line[3]     
        start2 = int(line[4])-1
        end2   = int(line[5])
        color = line[-1].split(",")[0].split("=")[-1]
        source = (name1, start1, end1, 500)
        destination = (name2, start2, end2, 500)
        circle.chord_plot(source, destination, edgecolor=color_map[color],linewidth=0.3)


circle.figure

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(phi, rho)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def offset(theta,r,wx,wy):
    x,y=pol2cart(r,theta)
    return cart2pol(x+wx,y+wy)


theta1= 0
r1=500+240*(0.87-vmin2)/(vmax2-vmin2)
theta1,r1 = offset(theta1,r1,-0,-85)
theta2= 0
r2=500+240*(1.87-vmin2)/(vmax2-vmin2)
theta2,r2 = offset(theta2,r2,-0,-85)
theta3= 0
r3=500+240*(2.87-vmin2)/(vmax2-vmin2)
theta3,r3 = offset(theta3,r3,0,-85)

plt.text(theta1,r1,"CN=1",rotation=0,fontsize=8)
plt.text(theta2,r2,"CN=2",rotation=0,fontsize=8)
plt.text(theta3,r3,"CN=3",rotation=0,fontsize=8)

plt.savefig(args.o)
