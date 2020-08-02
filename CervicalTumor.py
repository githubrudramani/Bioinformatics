#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 18:04:53 2020

@author: rudramani
"""
'''
Build machine-learning models from raw RNA-Seq count dat from
MLSeq package http://bioconductor.org/packages/release/bioc/html/MLSeq.html
we will work with the cervical count
data. Cervical data is from an experiment 
that measures the expression levels of 714 
miRNAs of human samples . There are 29 tumor
and 29 non-tumor cervical samples and these two 
groups can be treated as two separate classes for 
classification purpose. 
'''
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

'''
f = open('cervical.txt')
f.readline()
try:
    f = open("myfile", "r")
except IOError:
    print("the file myfile does not exists!!")
f.seek(0)
f.read()
'''

### Other Way the data is in dataframe structure
df = pd.read_csv('cervical.txt', sep ='\s+' )
df.shape
#(714, 58)
df.info()
df.describe()
df.isnull().sum()
df.head()
### converting samples as index and species as columns
df_Transpose = df.T
X = df_Transpose.values
### Labeling Tumor and non Tumor
m,n = df_Transpose.shape
y  =[0 for i in range(int(m/2))] + [1 for i in range(int(m/2))]



'''
Cervical data contains 714 miRNA mapped counts given in rows, 
belonging to 58 samples given in columns. First 29 columns 
of the data contain the miRNA mapped counts of non-tumor samples,
while the last 29 columns contain the count information of tumor samples.
'''
# Splitting the dataset into the Training set and Test set
from sklearn.cross_validation import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.25, random_state = 0)


# Feature Scaling
from sklearn.preprocessing import StandardScaler
sc = StandardScaler()
X_train = sc.fit_transform(X_train)
X_test = sc.transform(X_test)

"""
# Fitting SVM to the Training set
# Fitting SVM to the Training set
# Fitting SVM to the Training set
"""
from sklearn.svm import SVC
classifier = SVC(kernel = 'linear', random_state = 0)

classifier.fit(X_train, y_train)

# Predicting the Test set results
y_pred = classifier.predict(X_test)

# Making the Confusion Matrix
from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, y_pred)

from sklearn.metrics import accuracy_score
score = accuracy_score(y_test, y_pred)*100
print("Accuracy from SVM Model is: ", score)


"""
# Fitting Logistic Regression to the Training set
# Fitting Logistic Regression to the Training set
# Fitting Logistic Regression to the Training set
"""
from sklearn.linear_model import LogisticRegression
classifier = LogisticRegression(random_state = 0)
classifier.fit(X_train, y_train)

#Predictiong test set results
y_pred = classifier.predict(X_test)

# Making the Confusion matrix
from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, y_pred)
score2 = accuracy_score(y_test, y_pred)*100


"""
# Fitting  Random Forest classifier to the Training set
# Fitting  Random Forest classifier to the Training set
# Fitting  Random Forest classifier to the Training set
"""
from sklearn.ensemble import RandomForestClassifier
classifier = RandomForestClassifier(n_estimators = 10, criterion = 'entropy', \
                                    random_state =0)
classifier.fit(X_train, y_train)


# Predicting the Test set results
y_pred = classifier.predict(X_test)

# Making the Confusion Matrix
from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, y_pred)
score3 = accuracy_score(y_test, y_pred)*100
"""
# Fitting DT classifier to the Training set
# Fitting DT classifier to the Training set
# Fitting DT classifier to the Training set
"""
from sklearn.tree import DecisionTreeClassifier
classifier = DecisionTreeClassifier(criterion = 'entropy', random_state = 0)
classifier.fit(X_train, y_train)

# Predicting the Test set results
y_pred = classifier.predict(X_test)

# Making the Confusion Matrix
from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, y_pred)
score4 = accuracy_score(y_test, y_pred)*100

"""
# Fitting K-NN to the Training set
# Fitting K-NN to the Training set
# Fitting K-NN to the Training set
"""
from sklearn.neighbors import KNeighborsClassifier
classifier = KNeighborsClassifier(n_neighbors = 5, metric = 'minkowski', p = 2)
classifier.fit(X_train, y_train)


# Predicting the Test set results
y_pred = classifier.predict(X_test)

# Making the Confusion Matrix
from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, y_pred)
score5 = accuracy_score(y_test, y_pred)*100

"""
#
# Fitting Naive Bayes classifier to the Training set
#
"""
from sklearn.naive_bayes import GaussianNB
classifier = GaussianNB()
classifier.fit(X_train, y_train)

# Predicting the Test set results
y_pred = classifier.predict(X_test)

# Making the Confusion Matrix
from sklearn.metrics import confusion_matrix
cm = confusion_matrix(y_test, y_pred)

score6 = accuracy_score(y_test, y_pred)*100

"""
#
# Using Deep Learning Framework
#
"""

import keras
from keras.layers import Dense
from keras.models import Sequential
from keras.utils import to_categorical
import matplotlib.pyplot as plt
import tensorflow as tf 



y_train_cat = to_categorical(y_train)
y_test_cat = to_categorical(y_test)

sample, species = X_train.shape
num_classes = y_test_cat.shape[1]
print(num_classes)

def classification_model():
    # create mode
    model = Sequential()
    model.add(Dense(sample, activation = 'relu', input_shape = (species,)))
    model.add(Dense(25, activation = 'relu'))
    model.add(Dense(5, activation = 'relu'))
    #model.add(Dense(num_classes, activation = 'softmax')) ## need to upgrade
    model.add(Dense(2, activation = tf.nn.softmax))
    
    
    # compile model
    model.compile(optimizer = 'adam', loss = 'categorical_crossentropy',
                  metrics = ['accuracy'])
    return model


# build model
model = classification_model()

# fit the model
model.fit(X_train,y_train_cat, validation_split = 0.10, epochs = 10, verbose = 2)

# evaluate the model
scores = model.evaluate(X_test, y_test_cat, verbose = 0)
prediction = model.predict(X_test)
np.around(prediction, decimals = 0)

print('Accuracy: {}% \n Error: {}'.format(scores[1], 1 - scores[1]))  



print("Accuracy from SVM Model is: ", score)
print("Accuracy from logistic regression is: ", score2)
print("Accuracy from Decision Tree is: ", score4)
print("Accuracy from RandomForest is: ", score3)
print("Accuracy from KNN is: ", score5)
print("Accuracy from Naive Bayes  is: ", score6)
print("Accuracy from Neural Network  is: ", scores[1])



print("Hurray!!! The best Model in KNN")











