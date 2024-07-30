from pyjet import cluster
from pyjet.testdata import get_event
from pyjet import DTYPE_EP, DTYPE_PTEPM
from numpy.lib.recfunctions import append_fields
from numpy.testing import assert_array_equal
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


#Plotting
plt.style.use("./mystyle.mplstyle")
#Plots a step histogram to a fixed style
def StepHistoPlotter(dats, bins, xlim, labels, xlabel, title): 
    savedir="Plots/"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(dats)):
        dat, lab = dats[i], labels[i]
        plt.hist(dat, weights=np.zeros_like(dat)+1/dat.size, bins=bins, histtype='step',density=False,linewidth=1.5,label=lab)
    plt.xlabel(xlabel)
    plt.ylabel("Relative occurence")
    plt.xlim(xlim)
    handles, labels = ax.get_legend_handles_labels()
    new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]
    plt.legend(handles=new_handles, labels=labels)
    plt.savefig(savedir+title+".pdf", bbox_inches='tight')
    return

#Subjettiness ratios, just using the NSubjettines function
def NSubjetRatio(jet, R_subjet, N_numer, N_denom):
    return NSubjet(jet, R_subjet, N_numer)/NSubjet(jet, R_subjet, N_denom)

#Reads the given file and creates a list of jet constituents
def DatReader(filename):
    #Getting the full data
    datdir = "Data/"
    with open(datdir+filename+"_out.dat", 'r') as file:
        Dat = file.readlines()[3:]
    
    #Finding the start and end line for each jet, then extracting the data
    jets=[]
    jet_starts = [i for i in range(len(Dat)) if Dat[i][0] == "c"] #The jets are all prefaced by a line "constituents of jet n in event m" 
    N_jets = len(jet_starts)-1
    for i in range(N_jets):
        start, end = jet_starts[i]+1, jet_starts[i+1] #No +1 on end due to index slicing method below
        jet_parts = []
        
        #Getting the particle 4 momenta in the required form
        for line in Dat[start:end]:
            #Better way to do this that doesn't need the empty array - i.e. just go straight to the array?
            part  = np.empty(1, dtype=DTYPE_PTEPM)
            split = line.strip().split()
            part['pT'], part['eta'], part['phi'], part['mass'] = float(split[1]), float(split[2]), float(split[3]), float(split[4])
            jet_parts.append(part)
        
        #Appending to the jet list
        jets.append(np.array(jet_parts))
        if i % 10000 == 0: print("Jet",i+1,"read in")
        
    return jets
