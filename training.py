# -*- coding: utf-8 -*-

""" to check if tf is detecting your GPU un this python commands
import tensorflow as tf
physical_devices = tf.config.list_physical_devices()
print(physical_devices)
"""
#set file for species_redundancynumber or genus_redundancynumber
#options: genus_70, genus_80, genus_90, species_70, species_80, species_90
file = 'species_90'
filename = '../datasets/'+file+'.csv'

#set NDG to True to test with NDG dataset or False to test with 30% of the selected dataset
NDG = True

#model genus configurations
if file.split('_')[0] == 'genus':
    drop_rate = 0.35
    layer1 = 1000
    layer2 = 300
    hidden_activation = 'selu'
    output_activation = 'softmax'
    early_stop_delta=5.0e-3
    batchSize = 50

#model species configurations    
if file.split('_')[0] == 'species':
    drop_rate = 0.30
    layer1 = 3000
    layer2 = 200
    hidden_activation = 'selu'
    output_activation = 'sigmoid'
    early_stop_delta=5.0e-2
    batchSize = 30

log = 'results.log'
with open(log, 'w') as f:
    f.write('Layer 1: '+str(layer1)+'\n')
    f.write('Layer 2: '+str(layer2)+'\n')
    f.write('Output activation: '+output_activation+'\n')
    f.write('Early stopping: '+str(early_stop_delta)+'\n')
    f.write('NDG Test: '+str(NDG)+'\n')

from pandas import read_csv    
x_train = read_csv(filename, header = None)
x_train.drop(x_train.columns[-1], axis=1, inplace=True)
y_train = x_train.iloc[:, -1].tolist()
x_train.drop(x_train.columns[-1], axis=1, inplace=True)
x_train = x_train.values.tolist()
from sklearn.preprocessing import LabelBinarizer
LB = LabelBinarizer()
y_train = LB.fit_transform(y_train)
y_label = LB.classes_

if NDG == True:
    filename = file.split('_')[0]+'_NDG.csv'
    x_val = read_csv(filename, header = None)
    x_val.drop(x_val.columns[-1], axis=1, inplace=True)
    y_val = x_val.iloc[:, -1].tolist()
    x_val.drop(x_val.columns[-1], axis=1, inplace=True)
    x_val = x_val.values.tolist()
    y_val = LB.transform(y_val)
else:
    from sklearn.model_selection import train_test_split
    x_train, x_val, y_train, y_val = train_test_split(x_train, y_train, test_size=0.3)

from tensorflow.data import Dataset
train_dataset = Dataset.from_tensor_slices((x_train,y_train)).batch(batchSize,drop_remainder=True)
val_dataset = Dataset.from_tensor_slices((x_val,y_val)).batch(batchSize,drop_remainder=True)

from tensorflow import keras as k
model = k.Sequential()
model.add( k.layers.Dropout(rate = drop_rate))
model.add( k.layers.Dense(units = layer1, activation = hidden_activation, kernel_regularizer = k.regularizers.l2()))
model.add( k.layers.Dense(units = layer2, activation = hidden_activation))
model.add(k.layers.Dense(units=len(y_label),activation= output_activation))
model.compile( optimizer=k.optimizers.Adam(learning_rate=1e-4), loss = 'categorical_crossentropy', metrics=["accuracy"])

def ClassWeights( data ):
    from numpy import argmax, unique
    from sklearn.utils import class_weight
    data = argmax(data, axis=1)
    data = class_weight.compute_class_weight(class_weight = 'balanced', classes = unique(data), y = data)
    data = {i : data[i] for i in range(len(data))}
    return data

from tensorflow.keras.callbacks import EarlyStopping
early_stop = EarlyStopping(monitor='accuracy', patience=5, restore_best_weights=True, min_delta=early_stop_delta)

history = model.fit(train_dataset, validation_data = val_dataset, 
                    epochs=100, callbacks = [early_stop],
                    class_weight = ClassWeights(y_train), verbose=2,
                    use_multiprocessing=True)

y_pred = model.predict(x_val)
y_pred = LB.inverse_transform(model.predict(x_val))
y_val = LB.inverse_transform(y_val)
from sklearn.metrics import classification_report
metrics_report = classification_report(y_val,y_pred)
with open(log, '+a') as f:
    f.write('\n'+str(metrics_report))
    
def plot_confusion_matrix(dataTrue, dataPred, labels, normalize=True):
    from sklearn.metrics import confusion_matrix
    cm = confusion_matrix(dataTrue, dataPred, labels=labels)
    
    if normalize:
        cm = [[100*j/sum(i) for j in i if sum(i) != 0] for i in cm if None not in i]        
        cm = [i for i in cm if i] #Tirando listas vazias
    
    import matplotlib.pyplot as plt
    plt.figure(figsize=(23, 27))
    plt.imshow(cm, interpolation='nearest', cmap=plt.get_cmap('viridis'))
    plt.colorbar()
    tick_marks = numpy.arange(len(labels))
    plt.xticks(tick_marks, labels, rotation=90, fontsize=11)
    plt.yticks(tick_marks, labels, fontsize=11)
    plt.tight_layout()
    plt.ylabel('True label')
    plt.xlabel('Predicted label')
    plt.savefig('confusion_matrix.png', dpi=300)

  plot_confusion_matrix(y_val, y_pred, y_label)
