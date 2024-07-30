import os,sys
import json
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt
import networkx as nx
def read_json(filepath,max_events=0):
    all_events=[]
    with open(filepath,'r') as f:
        count=0
        for line in f:
            all_events.append(json.loads(line,encoding='ascii'))
            count+=1
            if count==max_events: break    
    #print (all_events,len(all_events))
    return all_events
    
   
all_vars=['p_m', 'p_pt', 'kt', 'psi', 's_pt', 'Delta', 'z', 'h_pt'] 



            
def get_lund_prim_graph(events,var_keys=None,max_nodes=15,log=None,return_sparse=True):
    '''get the feature vector and adjacency matrix of the primary Lund plane
            
            events: list of events each(jet) having a list of dictionaries with keys in all_vars,
                   each jet starts with the full jet(ind 0), next one(ind 1) is the harder 
                   pt daughter and the third one(ind 2) is the softer one,
                   from then on the odd indices correspond to the harder 
                   daughter of the previous odd integer and the even indices is 
                   softer daughter of the ind-3 particle
                
            returns: feature_matrix,adj_matrix,var_names
                   feature_matrix of shape (num_events,max_nodes,num_features), 
                   adj_matrix of shape (num_events,max_nodes,max_nodes)
                   var_names is the name of the feature variables
            return_sparse: if True converts the second and third axis into scipy.sparse.csr_matrix 
                           otherwise returns the numpy.ndarray
       currently returns adj_matrix per event to allow for future extension to full Lund Plane
       (use only one instance for the time being)
    '''
    assert max_nodes%2==1,'Max_nodes should be odd'
    if var_keys is None: var_keys=event[0,0].keys()
    if log is None: log=[False for _ in var_keys]
    adj_matrix=np.zeros((max_nodes,max_nodes),dtype='int32')
    #draw_graph(adj_matrix)
    ############ adjacency matrix ####################
    #Eadj_matrix[0]=np.array([
    for i in range(max_nodes): 
        if i==0:
            adj_matrix[0,0:3]=1
            continue          
        if i%2 !=0:
            adj_matrix[i,i-1:i+1]=1
        else:
            adj_matrix[i,i-2:i+2]=1
            adj_matrix[i,i-1]=0
    ##################################################
    #draw_graph(adj_matrix)
    '''
    adj_matrix=np.array([[1,1,1,0,0],
                         [1,1,0,0,0],
                         [1,0,1,1,0],
                         [0,0,1,1,0],
                         [0,0,1,0,1]])
                         '''
    feature_matrix=np.zeros((len(events),max_nodes,len(var_keys)),dtype='float32')
    points=np.zeros((len(events),max_nodes,2),dtype='float32')
    mask=np.ones((len(events),max_nodes,1))
    var_names=[]
    for key,l in zip(var_keys,log):
        if l:var_names.append('log_'+key)
        else:var_names.append(key)
    key_range=[_ for _ in range(len(var_keys))]
    all_events=[]
    num_consts=[]
    for event_index,jet in enumerate(events):
        num_consts.append(len(jet))
        mask[event_index,len(jet):]=0
        for i,node in enumerate(jet):
            #print (node)
            #adj_matrix[event_index,i,i:i+3]=1
            for ind,key,l in zip(key_range,var_keys,log):
                #print (key)
                if l:feature_matrix[event_index,i,ind]=np.log(node[key])
                else: feature_matrix[event_index,i,ind]=node[key]
        
        #print (feature_matrix[event_index],'\n',adj_matrix[event_index])
    if return_sparse:
        adj_matrix=np.array([sparse.csr_matrix(item) for item in adj_matrix])
        feature_matrix=np.array([sparse.csr_matrix(item) for item in feature_matrix])
    print ('Max number of constituents:',max(num_consts))
    #print (feature_matrix.shape,adj_matrix.shape)
    #adj_matrix=adj_matrix[num_consts.index(max(num_consts))]
    #print (adj_matrix)
    return feature_matrix,points,mask,adj_matrix,var_names
    
    
    
#return adj_train, adj, features_train, features, labels, idx_train, idx_val, idx_test
def draw_graph(adj_matrix,save_path=None):
    import networkx as nx
    print (adj_matrix)
    graph=nx.from_numpy_matrix(adj_matrix)
    nx.draw(graph)
    plt.show()
    
    
def load_data(filenames=None,**kwargs):
    for i,filename in enumerate(filenames):
        all_events=read_json(filename,max_events=kwargs.get('max_events',0))
        features,mask,adj_matrix,var_names=get_lund_prim_graph(all_events,
                                                        max_nodes=kwargs.get('max_nodes',31),
                                                        var_keys=kwargs.get('var_keys',['kt','Delta','z']),
                                                        log=kwargs.get('log',[True for _ in range(len(kwargs.get('var_keys',[0,0,0])))]),
                                                        return_sparse=kwargs.get('return_sparse',False))
        if i==0:
            allx=features
            allmask=mask
            ally=np.zeros((len(features),len(filenames)))
            ally[:,0]=1.
        else:
            allx=np.concatenate((allx,features),axis=0)
            allmask=np.concatenate((allmask,mask),axis=0)
            y=np.zeros((len(features),len(filenames)))
            y[:,i]=1.
            ally=np.concatenate((ally,y),axis=0)
    #draw_graph(adj_matrix)
    #sys.exit()
    np.random.seed(12)
    indices=np.arange(len(allx))
    np.random.shuffle(indices)
    train_split,test_end=int(0.6*len(allx)),int(0.2*len(allx))
    train_indices,test_indices,val_indices=indices[:train_split],indices[train_split:-test_end],indices[-test_end:]
    train_indices,test_indices,val_indices=np.sort(train_indices),np.sort(test_indices),np.sort(val_indices)
    x=allx[train_indices]
    mask=allmask[train_indices]
    y=ally[train_indices]
    if kwargs.get('model_name','particle_net')=='particle_net':
        val_x=allx[val_indices]
        val_y=ally[val_indices]
        val_mask=all_mask[val_indices]
        return 
    print (indices,train_indices,test_indices,val_indices)
    print (len(indices),len(train_indices),len(test_indices),len(val_indices))
    adj_matrix=nx.convert_matrix.from_numpy_array(adj_matrix)
    adj_matrix=nx.adjacency_matrix(adj_matrix)
    return adj_matrix,adj_matrix,np.matrix(x),np.matrix(allx),ally,train_indices,val_indices,test_indices
    
    
    draw_graph(adj_matrix)
    
    

            
if __name__=='__main__':load_data(filenames=['./jet_data/Fjets.json'],max_events=0)
 





