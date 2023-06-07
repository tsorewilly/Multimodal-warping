import tensorflow as tf
import numpy as np
from scipy import io
from numpy.linalg import norm
from itertools import product
import datetime
#from load_seq_data import load_train_data
from tensorflow.keras import Model
from tensorflow.keras import backend as K
from tensorflow.python.keras import backend as PK
from tensorflow.keras.optimizers import Adam, SGD
from tensorflow.keras.models import Sequential, load_model
from tensorflow.keras.layers import Dense, Dropout, BatchNormalization, Input, Conv1D, \
     UpSampling1D, Flatten, MaxPooling1D, GlobalAveragePooling1D, ReLU, Softmax
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.callbacks import ModelCheckpoint, ReduceLROnPlateau, EarlyStopping, TensorBoard
from losses import LossHistory, figPlot, trainPlot
import metrics as mtr
from tensorflow.keras.utils import plot_model
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from pycm import *

log_dir = "logs/"

from tensorflow.python.client import device_lib
print("Num GPUs Available: ", len(tf.config.experimental.list_physical_devices('GPU')))

#Load the preprocessed data from the matlab file
C = io.loadmat('warped_signals-6_Operators-17_15_trials.mat')
MM_Data = C['SpectSegMusc']

# extract dataset
Warped_Sequence_Data = MM_Data[:,0]
Warped_Instance_Data = MM_Data[:,1]
Warped_Labels = MM_Data[:,2]

#Set for network
Warped_Sequence = np.stack(Warped_Sequence_Data, axis=0)
Warped_Labels = to_categorical(Warped_Labels-1)


ModelName= 'BioCAS-CNN.v.1'
# CNN network design
def cnn_model():
    model = Sequential()
    input_shape = Input(shape = (Warped_Sequence.shape[1], Warped_Sequence.shape[2]), name='Input_Layer') #Sequence dimension: 100 x 26
    cnn_mod = Conv1D(filters=64, kernel_size=3, padding='same', activation='relu', name='Convolution_01')(input_shape)
    cnn_mod = BatchNormalization(name='Normalization_01')(cnn_mod)
    cnn_mod = MaxPooling1D(pool_size=4, strides=1, padding='valid', name='MaxPooling')(cnn_mod)
    cnn_mod = Dense(64, activation='softmax', name='Dense_01')(cnn_mod)
    cnn_mod = Dropout(0.45, name='DropOut_01')(cnn_mod)
    
    cnn_mod = Conv1D(filters=32, kernel_size=3, padding='same', activation='relu', name='Convolution_02')(cnn_mod)
    cnn_mod = BatchNormalization(name='Normalization_02')(cnn_mod)
    cnn_mod = Dense(64, activation='softmax', name='Dense_02')(cnn_mod)
    cnn_mod = MaxPooling1D(pool_size=4, strides=1, padding='same')(cnn_mod) #padding = 'valid' gave ValueError: Negative dimension size.
    cnn_mod = UpSampling1D(size=4, name='UpSampling_2')(cnn_mod)
    cnn_mod = Dropout(0.45, name='DropOut_02')(cnn_mod)
    
    output = Dense(6, activation='sigmoid', name='Dense_sigmoid')(cnn_mod)
    cnn_mod = Flatten(name='Flatten')(cnn_mod)
    cnn_mod = ReLU(negative_slope=0.01)(cnn_mod)
    
    output = Dense(2, activation='softmax', name='Dense_softmax')(cnn_mod) #activation='sigmoid'
    
    model = Model(inputs=[input_shape], outputs=[output], name='Output_Layer')
    plot_model(model, to_file='CNN-05-24.png', show_shapes=True, show_layer_names=True, dpi=100)
    
    # Save method, save once in 1 generation
    adam = Adam(lr=1e-4, decay=1e-6)
    model.compile(loss='mse', optimizer=adam, metrics=['acc', mtr.rmse, 'mae', 'mse'])#,  'cosine_proximity'])
    model.build(input_shape)
    
    return model


#define a function to fit the model
def fit_and_evaluate(t_x, val_x, t_y, val_y, EPOCHS=200, BATCH_SIZE=128):
    model = None
    model = cnn_model()
    strt_time = datetime.datetime.now()
    results = model.fit(t_x, t_y, epochs=EPOCHS, batch_size=BATCH_SIZE, 
                        shuffle=False, verbose=1, validation_split=0.1,
                        callbacks=[tbCallBack, checkpoint_period, reduce_lr, history, early_stopping])

    timedelta = datetime.datetime.now() - strt_time
    model_train_time = timedelta.total_seconds()
    print("Training completed in ", timedelta.total_seconds(), "s \n Saving the model.")
    
    plt.plot(results.epoch, np.array(results.history['val_loss']), label='Val loss')

    print("Val Score: ", model.evaluate(val_x, val_y))
    return results

checkpoint_period = ModelCheckpoint(log_dir + ModelName + '-model.h5',monitor='val_loss', save_weights_only=False, save_best_only=True, period=1)
reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.001, patience=3, verbose=1)
early_stopping = EarlyStopping(monitor='val_loss', min_delta=0, patience=10, verbose=1)
tbCallBack = TensorBoard(log_dir=log_dir+'./graph', histogram_freq=1, write_graph=True, write_images=True)
history = LossHistory()

Softmax_model = cnn_model()
Softmax_model.summary()

# Define the Keras TensorBoard callback.

n_folds = 5
epochs = 200
batch_size = 128
model_history = []

train_x, test_x, train_y, test_y = train_test_split(Warped_Sequence, 
                                                    Warped_Labels, 
                                                    random_state = 101, 
                                                   test_size=0.25)

for i in range(n_folds):
    print("Training on Fold: ",i+1)
    t_x, val_x, t_y, val_y = train_test_split(train_x, train_y, test_size=0.1, 
                                               random_state = np.random.randint(1,1000, 1)[0])
    
    model_history.append(fit_and_evaluate(t_x, val_x, t_y, val_y, epochs, batch_size))
    print("======="*12, end="\n\n\n")


plt.title('CNN: Accuracies vs Epochs')
for i in range(n_folds):
    plt.plot(model_history[i].history['acc'], label='Training Fold '+str(i))
    plt.plot(model_history[i].history['val_acc'], label='Training Fold '+str(i))
plt.legend()
plt.show()

# Final evaluation of the model
dependency = {'rmse': mtr.rmse}
#model = keras.models.load_model(self.output_directory + 'best_model.hdf5', custom_objects=dependencies)
model = load_model(log_dir+ModelName+'-model.h5', custom_objects=dependency)
model_prediction = model.predict(test_x)

y_pred = np.argmax(model_prediction, axis=1)
test_y = np.argmax(test_y, axis=1)

figPlot(test_y, y_pred, ax_fz=30, fz=30, saveAs = 'Output-CNN.png', title='CNN')


