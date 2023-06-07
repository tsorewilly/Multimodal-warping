import tensorflow.keras
from tensorflow.keras import backend as K


def loss(y_true, y_pred):
    return K.categorical_crossentropy(y_true, y_pred)

def rmse(y_true, y_pred):
    return K.sqrt(K.mean(K.square(y_pred - y_true), axis=-1))

# Write a Loss History class to save the loss and acc of the training set
# Of course, I can not do this at all, I can directly use the history object returned by the model.fit() method to do it
class LossHistory(tensorflow.keras.callbacks.Callback):
    def on_train_begin(self, logs={}):
        self.losses = {'batch': [], 'epoch': []}
        self.accuracy = {'batch': [], 'epoch': []}
        self.mse = {'batch': [], 'epoch': []}
        self.mae = {'batch': [], 'epoch': []}
        self.mape = {'batch': [], 'epoch': []}
        self.rmses = {'batch': [], 'epoch': []}
        self.cos_prox = {'batch': [], 'epoch': []}
        self.val_loss = {'batch': [], 'epoch': []}
        self.val_acc = {'batch': [], 'epoch': []}
        self.val_mse = {'batch': [], 'epoch': []}
        self.val_mae = {'batch': [], 'epoch': []}
        self.val_mape= {'batch': [], 'epoch': []}
        self.val_rmse = {'batch': [], 'epoch': []}
        self.val_cos_prox = {'batch': [], 'epoch': []}

    def on_batch_end(self, batch, logs={}):
        self.losses['batch'].append(logs.get('loss'))
        self.accuracy['batch'].append(logs.get('accuracy'))
        self.mse['batch'].append(logs.get('mean_squared_error'))
        self.mae['batch'].append(logs.get('mean_absolute_error'))
        self.mape['batch'].append(logs.get('mean_absolute_percentage_error'))
        self.rmses['batch'].append(logs.get('rmse'))
        self.cos_prox['batch'].append(logs.get('cosine_proximity'))
        self.val_loss['batch'].append(logs.get('val_loss'))
        self.val_acc['batch'].append(logs.get('val_accuracy'))
        self.val_mse['batch'].append(logs.get('val_mean_squared_error'))
        self.val_mae['batch'].append(logs.get('val_mean_absolute_error'))
        self.val_mape['batch'].append(logs.get('val_mean_absolute_percentage_error'))
        self.val_rmse['batch'].append(logs.get('val_rmse'))
        self.val_cos_prox['batch'].append(logs.get('val_cosine_proximity'))

    def on_epoch_end(self, batch, logs={}):
        self.losses['epoch'].append(logs.get('loss'))
        self.accuracy['epoch'].append(logs.get('accuracy'))
        self.mse['epoch'].append(logs.get('mean_squared_error'))
        self.mae['epoch'].append(logs.get('mean_absolute_error'))
        self.mape['epoch'].append(logs.get('mean_absolute_percentage_error'))
        self.rmses['epoch'].append(logs.get('rmse'))
        self.cos_prox['epoch'].append(logs.get('cosine_proximity'))
        self.val_loss['epoch'].append(logs.get('val_loss'))
        self.val_acc['epoch'].append(logs.get('val_accuracy'))
        self.val_mse['epoch'].append(logs.get('val_mean_squared_error'))
        self.val_mae['epoch'].append(logs.get('val_mean_absolute_error'))
        self.val_mape['epoch'].append(logs.get('val_mean_absolute_percentage_error'))
        self.val_rmse['epoch'].append(logs.get('val_rmse'))
        self.val_cos_prox['epoch'].append(logs.get('val_cosine_proximity'))
