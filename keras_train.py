#!/usr/bin/env python
# coding: utf-8

# In[1]:
import sys
import os
import numpy as np
import awkward


# In[2]:


import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] %(levelname)s: %(message)s')

from hep_ml.commonutils import map_on_cluster
from hep_ml.io.saver import Unpickle,Pickle
from hep_ml.genutils import check_dir

# In[3]:


def stack_arrays(a, keys, axis=-1):
    flat_arr = np.stack([a[k].flatten() for k in keys], axis=axis)
    return awkward.JaggedArray.fromcounts(a[keys[0]].counts, flat_arr)


# In[4]:


def pad_array(a, maxlen, value=0., dtype='float32'):
    x = (np.ones((len(a), maxlen)) * value).astype(dtype)
    for idx, s in enumerate(a):
        if not len(s):
            continue
        trunc = s[:maxlen].astype(dtype)
        x[idx, :len(trunc)] = trunc
    return x


# In[5]:


class Dataset(object):

    def __init__(self, filepath, feature_dict = {}, label='label', pad_len=100, data_format='channel_first'):
        self.filepath = filepath
        self.feature_dict = feature_dict
        if len(feature_dict)==0:
            feature_dict['points'] = ['part_etarel', 'part_phirel']
            feature_dict['features'] = ['part_pt_log', 'part_e_log', 'part_etarel', 'part_phirel']
            feature_dict['mask'] = ['part_pt_log']
        self.label = label
        self.pad_len = pad_len
        assert data_format in ('channel_first', 'channel_last')
        self.stack_axis = 1 if data_format=='channel_first' else -1
        self._values = {}
        self._label = None
        self._load()

    def _load(self):
        logging.info('Start loading file %s' % self.filepath)
        counts = None
        with awkward.load(self.filepath) as a:
            self._label = a[self.label]
            for k in self.feature_dict:
                cols = self.feature_dict[k]
                if not isinstance(cols, (list, tuple)):
                    cols = [cols]
                arrs = []
                for col in cols:
                    if counts is None:
                        counts = a[col].counts
                    else:
                        assert np.array_equal(counts, a[col].counts)
                    arrs.append(pad_array(a[col], self.pad_len))
                self._values[k] = np.stack(arrs, axis=self.stack_axis)
        logging.info('Finished loading file %s' % self.filepath)


    def __len__(self):
        return len(self._label)

    def __getitem__(self, key):
        if key==self.label:
            return self._label
        else:
            return self._values[key]
    
    @property
    def X(self):
        return self._values
    
    @property
    def y(self):
        return self._label

    def shuffle(self, seed=None):
        if seed is not None:
            np.random.seed(seed)
        shuffle_indices = np.arange(self.__len__())
        np.random.shuffle(shuffle_indices)
        for k in self._values:
            self._values[k] = self._values[k][shuffle_indices]
        self._label = self._label[shuffle_indices]
        
        
class LundData(object):
    def __init__(self,run_names=['dijet','w_fat_z_inv'],data_dir='./processed_events',key='lund_features',
                 ext='.pickle',**kwargs):
        self.run_names=run_names
        self.data_dir=os.path.abspath(data_dir)
        self.key=key
        self.ext=ext
        self.trainX,self.valX,self.testX=None,None,None
        self.trainY,self.valY,self.testY=None,None,None
        self.load_data()
    def load_data(self):
        for i,run_name in enumerate(self.run_names):
            events=Unpickle(run_name+self.ext,load_path=self.data_dir)
            if i==0:
                X=events[self.key]
                Y=np.zeros((len(X),len(self.run_names)))
                Y[:,i]=1.
            else:
                temp_X=events[self.key]
                temp_Y=np.zeros((len(temp_X),len(self.run_names)))
                temp_Y[:,i]=1.
                X=np.concatenate((X,temp_X),axis=0)
                Y=np.concatenate((Y,temp_Y),axis=0)
        np.random.seed(12)
        indices=np.arange(len(X))
        np.random.shuffle(indices)
        train_split,test_end=int(0.6*len(X)),int(0.2*len(X))
        train_indices,test_indices,val_indices=indices[:train_split],indices[train_split:-test_end],indices[-test_end:]
        train_indices,test_indices,val_indices=np.sort(train_indices),np.sort(test_indices),np.sort(val_indices)
        self.trainX,self.valX,self.testX=X[train_indices],X[val_indices],X[test_indices]
        self.trainY,self.valY,self.testY=Y[train_indices],Y[val_indices],Y[test_indices]
        
        


# In[6]:


#train_dataset = Dataset('converted/test_file_0.awkd', data_format='channel_last')
#val_dataset = train_dataset#Dataset('converted/val_file_0.awkd', data_format='channel_last')

data=LundData()
# In[ ]:





# In[7]:


import tensorflow as tf
from tensorflow import keras
from tf_keras_model import get_particle_net, get_particle_net_lite


# In[8]:


model_type = 'particle_net_lite' # choose between 'particle_net' and 'particle_net_lite'
num_classes = 2#train_dataset.y.shape[1]
input_shapes ={'points':(30,3)}#train_dataset.X}#{'points':(30,2),'features':(30,3),'mask':(30,1)} #
print (input_shapes)
#input('press enter to continue:')
if 'lite' in model_type:
    model = get_particle_net_lite(num_classes, input_shapes)
else:
    model = get_particle_net(num_classes, input_shapes)
model.summary()
# In[9]:


# Training parameters
batch_size = 100 if 'lite' in model_type else 384
epochs = 30


# In[10]:


def lr_schedule(epoch):
    lr = 1e-3
    if epoch > 10:
        lr *= 0.1
    elif epoch > 20:
        lr *= 0.01
    logging.info('Learning rate: %f'%lr)
    return lr


# In[11]:


model.compile(loss='categorical_crossentropy',
              optimizer='nadam',#keras.optimizers.Adam(learning_rate=lr_schedule(0)),
              metrics=['accuracy'])
model.summary()


# In[12]:


# Prepare model model saving directory.
import os
save_dir = 'model_checkpoints'
model_name = '%s_model.{epoch:03d}.hdf5' % model_type
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)
filepath = os.path.join(save_dir, model_name)

# Prepare callbacks for model saving and for learning rate adjustment.
checkpoint = keras.callbacks.ModelCheckpoint(filepath=filepath,
                             monitor='val_acc',
                             verbose=1,
                             save_best_only=True)

lr_scheduler = keras.callbacks.LearningRateScheduler(lr_schedule)
progress_bar = keras.callbacks.ProgbarLogger()
tboard=keras.callbacks.TensorBoard(
    log_dir='./logs', histogram_freq=300, write_graph=True,
    write_images=True, update_freq=100, profile_batch=2,
    embeddings_freq=0, embeddings_metadata=None
)
callbacks = [tboard]
# In[13]:


#train_dataset.shuffle()
history=model.fit(data.trainX,data.trainY,
          batch_size=batch_size,
#           epochs=epochs,
          epochs=50, # --- train only for 1 epoch here for demonstration ---
          validation_data=(data.valX, data.valY),
          shuffle=True,
          callbacks=callbacks)
          
history=history.history
hist_path=check_dir('./training_history')
name='history_train'+str(len(os.listdir(hist_path)))
Pickle(history,name,save_path=hist_path)


# In[ ]:




