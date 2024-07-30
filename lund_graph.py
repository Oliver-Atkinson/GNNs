import os
import sys
import pickle

from pyjet import cluster
from pyjet.testdata import get_event
import numpy as np
from math import sqrt,log,atan
from ROOT import TLorentzVector

#from hep_ml.genutils import print_events,merge_flat_dict,pool_splitter,check_dir
#from hep_ml.io.saver import Unpickle,Pickle
#from hep_ml.hep.utils.fatjet import FatJet 

'''
 f=FatJet(tower[:,0:3],algorithm="antikt",r=1.,pt_min=200.)
        fatjets=f.Get()
        if len(fatjets)<1: continue
        cut_flow["fatjet"]+=1
        fatjet=fatjets[0]
        if abs(fatjet.eta)>2.5: continue
        cut_flow["eta"]+=1
        #if fatjet.pt <160. : continue
        #print (fatjet)
        #cut_flow["pt"]+=1
'''


AlgorithmDict={"antikt":-1,"CA":0,"kt":1}
def construct_vector(jet,maximum_particles=None):
    jet=np.array([[item.pt,item.eta,item.phi,item.mass] for item in jet])
    vectors=[]
    for item in jet[:maximum_particles]:
        vectors.append(np.array((item[0],item[1],item[2],item[3]),dtype=[('pT', 'f8'), ('eta', 'f8'), ('phi', 'f8'), ('mass', 'f8')]))
    return np.array(vectors)

#print (construct_vector(jet),jet.constituents())
#length=len(jet.constituents())
#for i,const in enumerate(jet.constituents()):
#    print (type(const),dir(const))
    #const.set_user_index(i)
    #break
#sys.exit()
seq=[]
hard=[]
soft=[]
constituents=[]

def lund_vars(hard,soft,log_keys=['kt','z','m_ab','delta'],jet_vec=None):
    h_vec,s_vec=TLorentzVector(hard.px,hard.py,hard.pz,hard.e),TLorentzVector(soft.px,soft.py,soft.pz,soft.e)
    return_dict={}
    return_dict['delta']=h_vec.DeltaR(s_vec)
    return_dict['kt']=s_vec.Pt()*return_dict['delta']
    return_dict['z']=soft.pt/(hard.pt+soft.pt)
    return_dict['m_ab']=(h_vec+s_vec).M()
    #return_dict['psi']=atan((h_vec.Eta()-s_vec.Eta())/(h_vec.Phi()-s_vec.Phi()))
    for key in log_keys: 
        #print (key)
        return_dict[key]=log(return_dict[key])
    if jet_vec is not None:
        return_dict['eta']=jet_vec.Eta()-s_vec.Phi()
        return_dict['phi']=jet_vec.DeltaPhi(s_vec)
    return return_dict
    

all_seq=[]
def get_sequence(start,coordinate_funcs=['kt','z','delta'],pad=30):
    seq=[]
    hard=[]
    soft=[]
    lund_coordinates=[]
    eta_phi=[]
    count=1
    jet=TLorentzVector(start.px,start.py,start.pz,start.e)
    while start.parents is not None:        
        assert start.parents[0].pt > start.parents[1].pt
        try: _=lund_vars(start.parents[0],start.parents[1],jet_vec=jet)
        except Exception as e:
            print (e,'Ignoring!')
            continue
        seq.append(start.parents)
        eta_phi.append([_.pop('eta'),_.pop('phi')])
        lund_coordinates.append([_[key] for key in coordinate_funcs])
        hard.append(start.parents[0])
        #print (count,' ',start,start.parents)
        soft.append(start.parents[1])
        start=start.parents[0]
        #if start.parents is None: const=start
        count+=1 
    #print (len(seq),len(hard),len(soft))
    assert len(lund_coordinates)<=pad
    while len(lund_coordinates)<pad: lund_coordinates.append([0 for _ in coordinate_funcs])
    return {'seq':seq,'hard':hard,'soft':soft,'lund_features':np.array(lund_coordinates),'points':np.array(eta_phi)}
    

def get_lund_coordinates(events):
     #print_events(events)
     f=FatJet(events['Tower'][0][:,0:3],algorithm="antikt",r=1.,pt_min=200.)
     lund_features=[]
     points=[]
     for i,fatjet in enumerate(events['FatJet']):
         cluster_seq,length=f.Recluster(fatjet,r=1.,algorithm='CA',return_val='sequence')
         #print (cluster_seq)
         start=cluster_seq.exclusive_jets(1)[0]
         temp_dict=get_sequence(start)
         lund_features.append(temp_dict['lund_features'])
         points.append(temp_dict['points'])
     lund_features=np.array(lund_features)
     points=np.array(points)
     return {'lund_features':lund_features,'points':points}
 
"""
store_path=check_dir('./processed_events')   
for run_name in ['short_dijet','dijet','w_fat_z_inv'][1:]:
    events=Unpickle(run_name+'.pickle',load_path='./passed_events')
    #for key,val in events.items():
    #    if len(val)>125000: events[key]=val[:125000]
    lund_coordinates=pool_splitter(get_lund_coordinates,events)
    print_events(lund_coordinates,name=run_name)
    #break
    Pickle(lund_coordinates,run_name+'.pickle',save_path=store_path)
"""



