import matplotlib.pyplot as plt
import os
import sys
import numpy as np
from numpy import linalg as LA
import skfda
from skfda.preprocessing.dim_reduction import FPCA

directoryName="Data"

def Average_function_substraction(ma):
    Aveg=[]
    for i in range(len(ma[0])):
        summation=[]
        for f in ma:
            summation.append(f[i])
        ave=sum(summation) / len(summation)
        Aveg.append(ave)
        for f in ma:
            f[i]=f[i]-ave
    return Aveg, ma
        
def Cov_matrix(ma_data):
    Cov_matrix=[]
    ma_data=ma_data.T
    for x1 in ma_data:
        cov_column=[]
        for y2 in ma_data:
            print(len(y2))
            cov_column.append(x1.dot(y2)/len(y2))
        Cov_matrix.append(cov_column)
    return np.array(Cov_matrix)
            
def FactorVectors(Energy_list, Matrix_data):
    structure_val_vector = []
    avg, Matrix_avg_data=Average_function_substraction(Matrix_data)
    Matrix_cov=Cov_matrix(Matrix_avg_data)
    #print(Matrix_cov)
    eig_values, eig_vectors = LA.eig(Matrix_cov)
    Factor_eigValues = eig_values/sum(eig_values)
    for i in range(len(Factor_eigValues)):
        structure_val_vector.append([Factor_eigValues[i],eig_vectors[i]])
    return Factor_eigValues, structure_val_vector, avg

def InnerProduct(x_axis,f1,f2):
    result=0.0
    for i in range(len(x_axis)):
        if i==0:
            result+=0.008*f1[i]*f2[i]
        result+=abs(x_axis[i]-x_axis[i-1])*f1[i]*f2[i]
    return result
    
def Clipping_CollectingData(FileName):
    file1 = open(directoryName+"/"+FileName,'r')
    data=[]
    for line in file1:
        try:
            float(line.split('\t')[0])
        except ValueError:
           continue
        dataList=line.split('\n')[0].split('\t')
        #print(dataList)
        if len(data)==0:
            for i in range(len(dataList)):
                data.append([])
        for j in range(len(data)):
            data[j].append(float(dataList[j]))
            #print(dataList[j])
    file1.close()
    MatrixData=data[1:]
    return data[0],MatrixData
    
ifExistDirectory = os.path.exists("./"+directoryName)
if not(ifExistDirectory):
    print("Directory "+directoryName+" does not exist")
    os.mkdir("./"+directoryName)
    print("Directory "+directoryName+" has been created")
    print("Please move the cycling text files you want to work with in this directory")
    sys.exit()
  
files=os.listdir("./"+directoryName)
if len(files)==0:
    print(directoryName+" is empty")
    print("Please move the cycling text files you want to work with in this directory")
    sys.exit()
    
for x in files:
    if x[-3:]!='txt':
        next
    else:
        Energy_list, MatrixD = Clipping_CollectingData(x)
        fd = skfda.FDataGrid(
            data_matrix=MatrixD,
            grid_points=Energy_list,
            )
        Avg_num,_ = Average_function_substraction(MatrixD)
        print(Avg_num)
        fpca_discretized = FPCA(n_components=5)
        fpca_discretized.fit(fd)
        scores=fpca_discretized.fit_transform(fd)
        components=fpca_discretized.components_
        avg_fun=fpca_discretized.mean_.data_matrix.reshape(len(Energy_list),)
        print(avg_fun)
        modi=avg_fun
        for index in range(len(scores[0])):
            modi+=scores[0][index]*components.data_matrix[index].reshape(len(Energy_list),)
        #print(components0.data_matrix[0].reshape(len(Energy_list),))
        #components=fpca_discretized.singular_values_
        #print(components)
        #factors_value=fpca_discretized.explained_variance_ratio_
        #for i in components.data_matrix[0]:
        #    print(i[0])
        #f_final=Energy_list
        #for j in range(len(coeff)):
        #    tota_factors+=abs(factors_value[j])
        #    f_final+=coeff[j]*factors_vectors[j][1].real
        #print('Using {}% of variance'.format(tota_factors*100))
        plt.plot(Energy_list, Avg_num, label='Avg. num')
        plt.plot(Energy_list, avg_fun, label='Avg. library')
        plt.plot(Energy_list, Avg_num-avg_fun, label='Avg. diff')
        #plt.plot(Energy_list, MatrixD[0], label='Data[0]')
        #plt.plot(Energy_list, Avg_function+f_final, label='Aprox')
        plt.legend(loc='best')
        plt.show()
        

