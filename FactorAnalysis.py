import matplotlib.pyplot as plt
import os
import sys
import numpy as np
from numpy import linalg as LA
import math

directoryName="Data"

def saveData(x_axis,f1,f2,fileName,f_var):
    FileName_toSave="./"+directoryName+"/"+fileName+"_FPCA.txt"
    with open(FileName_toSave,"a") as txt_file:
        txt_file.write("Total Var. FPCA: "+str(f_var)+"\n")
        txt_file.write("Energy Loss\tNum. Avg.\tFPCA Avg.\n")
        for i in range(len(x_axis)):
            txt_file.write("\t".join([str(x_axis[i]),str(f1[i]),str(f2[i])]))
            txt_file.write("\n")

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
            cov_column.append(x1.dot(y2)/len(y2))
        Cov_matrix.append(cov_column)
    return np.array(Cov_matrix)
            
def FactorVectors(Energy_list, Matrix_data):
    structure_val_vector = []
    avg, Matrix_avg_data=Average_function_substraction(Matrix_data)
    Matrix_cov=Cov_matrix(Matrix_avg_data)
    eig_values, eig_vectors = LA.eig(Matrix_cov)
    Total_norm=0
    for i in range(len(eig_values)):
        eig_vectors[i].real=eig_vectors[i].real/math.sqrt(InnerProduct(Energy_list,eig_vectors[i].real,eig_vectors[i].real))
        structure_val_vector.append([eig_values[i],eig_vectors[i]])
        Total_norm+=abs(eig_values[i])
    Factor_eigValues=eig_values/Total_norm
    return Factor_eigValues, structure_val_vector, avg

def InnerProduct(x_axis,f1,f2):
    result=0.0
    for i in range(len(x_axis)):
        if i==0:
            result+=0.008*f1[i]*f2[i]
        result+=abs(x_axis[i]-x_axis[i-1])*f1[i]*f2[i]
    return result

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

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
    MatrixData=np.array(data[1:])
    #print(MatrixData.shape)
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
        for f in MatrixD:
            f=smooth(f,5)
        factors_total, factors_val_vectors, Avg_function = FactorVectors(Energy_list, MatrixD)
        _, Matrix_plot = Clipping_CollectingData(x)
        tota_factors=0
        for i in range(len(factors_total)):
            print('factor {}: {}'.format(i+1,abs(factors_total[i])))
            if abs(factors_total[i])<1.0*pow(10,-16):
                NUm_factors=i
                break
        print('Enter numbers of factors to build final data_scan:')
        Num = int(input())
        coeff=[]
        for f in Matrix_plot:
            co=[]
            g=f-Avg_function
            for j in range(Num):
                co.append(InnerProduct(Energy_list,g,factors_val_vectors[j][1].real))
            coeff.append(co)
        coeff=np.array(coeff)
        for j in range(Num):
            tota_factors+=abs(factors_total[j])
        #print('Using {}% of variance'.format(tota_factors*100))
        final_fun = 0.0*MatrixD[0] #gives the right dimension to the final object
        scan_num=2
        for n in range(len(coeff[scan_num-1])):
            final_fun+=coeff[scan_num-1][n]*factors_val_vectors[n][1].real
        #mag=InnerProduct(Energy_list,Matrix_plot[scan_num-1]-Avg_function,Matrix_plot[scan_num-1]-Avg_function)
        final_fun1=Matrix_plot[scan_num-1]-final_fun
        saveData(Energy_list,Avg_function,final_fun1,x[:-4],tota_factors)
        ##Plotting fig1
        fig1 = plt.figure()
        axes1=fig1.add_subplot(1,1,1)
        axes1.set_xlabel('Score fPCA1')
        axes1.set_ylabel('Score fPCA2')
        axes1.scatter(coeff.T[0], coeff.T[1], label='Scores', marker='o')
        axes1.legend(loc ='best')
        fig1.savefig("./"+directoryName+"/"+x[:-4]+'_scores'+'.png', dpi=300, facecolor = 'w')
        ##Plotting fig2
        fig2 = plt.figure()
        axes2=fig2.add_subplot(1,1,1)
        axes2.set_xlabel('Energy Loss (eV)')
        axes2.set_ylabel('Intensity (arb. units)')
        axes2.plot(Energy_list, smooth(final_fun1,5), label='FPCA Avg. using {}% var.'.format(tota_factors*100))
        axes2.plot(Energy_list, Avg_function, label='Numerical Avg.')
        axes2.legend(loc ='best')
        fig2.savefig("./"+directoryName+"/"+x[:-4]+'.png', dpi=300, facecolor = 'w')
        #plt.show()
        

