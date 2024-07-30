from Functions import *

#Desired output is a set of substructure variables
#For now focus on N-subjettiness - see 1011.2268 for a definition

#Takes a jets constituents and returns the N-subjettiness
def NSubjettiness(parts, R_subjet, N):
    #Clustering Particle 4-vectors of the form (pt, eta, phi, m) to N subjets
    assert len(parts) >= N, 'Too few particles in the jet'
    #anti-kt (p=-1) is not implemented in fastjet and only returns the N lowest pt jets
    #Include jet axes jiggling for proper axes minimisation as in fjcontrib (see notes)?
    subjets = cluster(parts,R=R_subjet,p=1,algo='genkt').exclusive_jets(N)
    
    #Initialising sum parts
    d0, tau_sum = 0, 0
    #Iterating over each particle in the full jet
    for part in parts:
        d0 += part['pT']
        DRs = [np.sqrt((part['eta']-subj.eta)**2 + (part['phi']-subj.phi)**2) for subj in subjets]
        tau_sum += part['pT']*min(DRs)
    d0 *= 1.5 #Original jet radius - update to match proper value (1.2?) - INPUT REQUIRED
    tau = tau_sum/d0
    return tau

#Taking a list of jets and clustering to microjets
def MicroJets(jets,R_micro,p,ptmin):
    #What about unclustered particles in the subjettines calculation?
    return [np.array(cluster(jet,R=R_micro,p=p,algo='genkt').inclusive_jets(ptmin=ptmin)) for jet in jets]

#Calculating the N subjettiness' of a list of jets
#Maybe put a loop over N in here
def SubjettinessCalc(jets,R_subj):
    tau1s, tau2s, tau3s = [], [], []
    #Going over each jet, getting the microjet information and using these to calculate the N subjettiness
    for jet in jets:
        constits = []
        for microjet in jet:
            part = np.empty(1, dtype=DTYPE_PTEPM) #Better way to do this that doesn't need the empty array - i.e. just go straight to the array?
            part['pT'], part['eta'], part['phi'], part['mass'] = microjet.pt, microjet.eta, microjet.phi, microjet.mass
            constits.append(part[0])
        
        if len(constits) >= 3: #Think more about this condition - Input required
            tau1s.append(NSubjettiness(np.array(constits), R_subj, 1))
            tau2s.append(NSubjettiness(np.array(constits), R_subj, 2))
            tau3s.append(NSubjettiness(np.array(constits), R_subj, 3))
    tau1s, tau2s, tau3s = np.array(tau1s), np.array(tau2s), np.array(tau3s)
    #Checks - maybe implement for the ratios too
    Nbigtau1s, Nbigtau2s, Nbigtau3s = sum(np.where(tau1s > 1, 1, 0)), sum(np.where(tau2s > 1, 1, 0)), sum(np.where(tau3s > 1, 1, 0))
    if Nbigtau1s != 0: print("Warning: "+'{:.3g}'.format(100*(Nbigtau1s/len(tau1s)))+"% of the jets have tau_1 > 1")
    if Nbigtau2s != 0: print("Warning: "+'{:.3g}'.format(100*(Nbigtau2s/len(tau2s)))+"% of the jets have tau_2 > 1")
    if Nbigtau3s != 0: print("Warning: "+'{:.3g}'.format(100*(Nbigtau3s/len(tau3s)))+"% of the jets have tau_3 > 1")

    return tau1s, tau2s, tau3s

#Actually doing the calculations
def Calculator(DatName,save=0,read=1, R_micro=0.1,p_micro=-1,ptmin_micro=5, R_subj=0.6): #INPUT REQUIRED ON R_subj
    print("Now performing subjet calculations for the "+DatName+" data")
    #Getting the full jet constituents, either from file or by calculating
    if read == 1:
        with open("Data/"+DatName+'_jets.pkl', 'rb') as f:
            jets = pickle.load(f)
        print("Jets read in")
    if save == 1: 
        jets = DatReader(DatName)
        with open("Data/"+DatName+'_jets.pkl', 'wb') as f:
            pickle.dump(jets, f)        
    
    #Clustering to microjets
    microjets = MicroJets(jets,R_micro,p_micro,ptmin_micro)
    n_micros = np.array([len(micros) for micros in microjets])
    print("Microjets found")
    
    #N-Subjettiness
    tau1s, tau2s, tau3s = SubjettinessCalc(microjets, R_subj)
    tau21s, tau32s = tau2s/tau1s, tau3s/tau2s
    
    return np.array([jets, microjets, n_micros, tau1s, tau2s, tau3s, tau21s, tau32s])


#Clear issues with some of the taus being > 1 - this is even more pronounced for the ratios

#A class would be far better here - implement this
tau_dict = {"QCD": Calculator("QCD"),
           "WZ": Calculator("WZ"),
           "Tops": Calculator("Tops"),
           "BSM": Calculator("BSM")}
id_dict = {3:[0.5,"1"],4:[0.1,"2"],5:[0.005,"3"],6:[1,"21"],7:[1,"32"]} #This is disgusting, resolve it
labels = ["QCD", r'$W^\pm$', r'$t$', r'$\phi$']

N_zonebins = 50
for i in range(2,8): #Update for full substructure stuff
    dat = [tau_dict[name][i] for name in tau_dict.keys()]
    if i == 2:
        StepHistoPlotter(dat, np.arange(0.5,21.5,1), [0.5,20], labels, "Microjets", "Nmicros")
    else:
        bmax = id_dict[i][0]
        Nbins = N_zonebins/bmax
        bins = np.arange(0,1+1/Nbins,1/Nbins)
        StepHistoPlotter(dat, bins, [0,bmax], labels, r'$\tau_{'+id_dict[i][1]+'}$', "Micros_tau"+id_dict[i][1])
    