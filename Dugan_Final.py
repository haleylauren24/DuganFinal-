#This code analyzes mass spectrometry data generated using Skyline,
#a program that extracts and organizes raw data taken directly from
#the mass spectrometer. The csv file provided contains mass spec data
#for proteins isolated from the liver of a Sirt5 knock-out mouse.
#Sirt5 is an enzyme that acts as a deacetylase, desuccinylase, and demalonylase
#to make post-translational modifications on lysine residues of proteins.
#Sirt5 has been shown to be important for mitochondrial metabolism
#and is relevant for understanding metabolic diseases.
#The csv file contains the fragmentation data for WT and Sirt5 KO proteins
#and shows the differences in post-translational modifications in the
#modifed sequence column.
#In this code, an interface will prompt the user to enter the number of 
#custom experimental conditions (2), WT and KO, and the code creates a tabular
#statistical report that includes the average peak area, standard deviation,
#and coefficient of variation.


import numpy as np
from scipy import stats
import scipy.stats as st
import sys
import csv
import os



fn = "Dugan_MSdata.csv"
    
with open(fn, 'r') as myfile:
    skylineFile = csv.reader(myfile, dialect='excel')
    data = []
    for row in skylineFile:
        data.append(row)
    myfile.close()
    

#Defining all of the necessary information required to create the output file 
#The output file will calculate the mean, sd, and cv for the specified 
#experimental conditions KO and WT based on the data file provided 

AreaCol = data[0].index("Area")
Areas2 = [row[AreaCol] for row in data[1 :]]
x=0
Areas=[]
while x<len(Areas2):
    if Areas2[x] == '#N/A':
        Areas.append('N/A')

    else:
        Areas.append(float(Areas2[x]))
    x+=1

#Defining the column with the fragment ion data 
ionsCol = data[0].index("FragmentIon")
ions = [row[ionsCol] for row in data[1 :]]

#Defining the column with the peptide names 
peptidesCol = data[0].index("ModifiedSequence")
peptides = [row[peptidesCol] for row in data[1 :]]

#Defining the column with the precursor M/Z ratio 
lightprecMZCol = data[0].index("PrecursorMz")
lprecs = [row[lightprecMZCol] for row in data[1 :]]

#Defining the column with the isotope distance rank 
idrCol = data[0].index("IsotopeDistRank")
idr = [row[idrCol] for row in data[1 :]]

#Defining the column with the replicate name 
rname = data[0].index("ReplicateName")
repname = [row[rname] for row in data[1 :]]

#Defining the column with the modified seq in the KO condition 
modseq = data[0].index("ModifiedSequence")
mseq = [row[modseq] for row in data[1 :]]

#Defining the column with the charge Z of the precursor 
prch = data[0].index("Z")
precch = [row[prch] for row in data[1 :]]

#Defining the column with the protein name 
protname = data[0].index("ProteinName")
proname = [row[protname] for row in data[1 :]]

#Defining the retention time
rettime = data[0].index("Avg RT")
rtime = [row[rettime] for row in data[1 :]]

#Defining the column with the description for the protein family
protDesc = data[0].index("ProteinDescription")
pdes = [row[protDesc] for row in data[1 :]]



i=0
j=0
k=1
z=0
KO=[]
WT=[]
posList=[]
indexList=[]
precAr=[[]]
repNameAr=[[]]
pepAr=[[]]
sucAr=[[]]
ionAr=[[]]
noteAr=[[]]
rankAr=[[]]
mseqAr=[[]]
precMzAr=[[]]
proNameAr=[[]]
precChAr=[[]]
rtAr=[[]]
prodesAr=[[]]

conditions = int(input('How many experimental conditions are there? '))
x=0
condAr=[]
while x<conditions:
    conds =input('Enter condition ' + str(x+1) +': ' )
    condAr.append(conds)
    x+=1

#Splitting up the peak areas and corresponding parameters from Skyline
#into a list of lists based on the product mass to charge
#ratio because this is what varies between the transitions.
    
while k<len(lprecs):
    if all([lprecs[j] == lprecs[k], mseq[j] == mseq[k], precch[j] == precch[k]]):
        precAr[z].append(Areas[j])
        repNameAr[z].append(repname[j])
        pepAr[z].append(peptides[j])
        ionAr[z].append(ions[j])
        mseqAr[z].append(mseq[j])
        rankAr[z].append(idr[j])
        precMzAr[z].append(lprecs[j])
        proNameAr[z].append(proname[j])
        precChAr[z].append(precch[j])
        rtAr[z].append(rtime[j])
        prodesAr[z].append(pdes[j])
    elif any([lprecs[j] != lprecs[k], mseq[j] != mseq[k], precch[j] != precch[k]]):
        precAr[z].append(Areas[k-1])
        repNameAr[z].append(repname[k-1])
        pepAr[z].append(peptides[k-1])
        ionAr[z].append(ions[k-1])
        mseqAr[z].append(mseq[k-1])
        rankAr[z].append(idr[k-1])
        precMzAr[z].append(lprecs[k-1])
        proNameAr[z].append(proname[k-1])
        precChAr[z].append(precch[k-1])
        rtAr[z].append(rtime[k-1])
        prodesAr[z].append(pdes[k-1])
        precAr.append([])
        repNameAr.append([])
        pepAr.append([])
        ionAr.append([])
        mseqAr.append([])
        rankAr.append([])
        precMzAr.append([])
        proNameAr.append([])
        precChAr.append([])
        rtAr.append([])
        prodesAr.append([])
        z+=1
    else:
        break
    j+=1
    k+=1

while k<len(lprecs):
    l=0
    if lprecs[j] == lprecs[k]:
        precAr[z].append(Areas[j])
        if k == len(lprecs)-1:
            precAr[z].append(Areas[k])
    elif lprecs[j] != lprecs[k]:
        precAr[z].append(Areas[k-1])
        precAr.append([])
        z+=1
    elif k == len(lprecs)-1:
        precAr[z].append(Areas[k])
        z+=1
    else:
        break
    j+=1
    k+=1

#p-value function
def t_test(array1, array2):
    num1 = len(array1)
    var1 = np.var(array1, ddof=1)
    num2 = len(array2)
    var2 = np.var(array2, ddof=1)
    degf = ((var1/num1 + var2/num2)**(2.0))/((var1/num1)**(2.0)/(num1-1) + (var2/num2)**(2.0)/(num2-1))
    t = (np.mean(array1) - np.mean(array2)) / np.sqrt(var1/num1 + var2/num2)

    pvalue = 1.0 - ( st.t.cdf(np.abs(t),degf) - st.t.cdf(-np.abs(t),degf) )    

    return pvalue





i=0
WT=[]
KO=[]
KOthenWT=[]
assign=[]
k=0
assign2=[]
i=0

#Based on the description provided, the information is split up into a list
#of lists with the KO data first and WT following.
while i<len(precAr):
    assign2=[]
    k=0
    while k<len(condAr):
        j=0
        assign=[]
        while j<len((precAr[i])):
            
            if condAr[k] in repNameAr[i][j]:
                assign.append(precAr[i][j])  
              
            
            j+=1
            
        assign2.append(assign)
        k+=1
    KOthenWT.append(assign2)
    i+=1

w=0
KOthenWT2=[]
while w<len(KOthenWT):
    y=0
    kos=[]
    while y<len(KOthenWT[w]):
        x=0
        kos2=[]
        while x<len(KOthenWT[w][y]):
            #if the area in the input file is not equal to N/A or 0 then append it
            if all ([KOthenWT[w][y][x] != 'N/A', KOthenWT[w][y][x] != 0]):
                kos2.append(KOthenWT[w][y][x])
            x+=1
        kos.append(kos2)
        y+=1
    KOthenWT2.append(kos)
    w+=1

MnWT=[]
SdWT=[]
CvWT=[]
MnKO=[]
SdKO=[]
CvKO=[]
w=0
k=0
ratioAr=[]
pvalueAr=[]

#calculating the mean, sd, and cv for the data of KO and WT
while w<len(KOthenWT2):
    y=0
    while y<len(KOthenWT2[w]):
        for num in KOthenWT2[w][y]:
            if num == 0:
                KOthenWT2[w][y].remove(num)
        wtmean = np.mean(KOthenWT2[w][y])
        MnWT.append(wtmean)
        sdArray=[]
        r=0
        while r <len(KOthenWT2[w][y]):
            sdArray.append(np.square((KOthenWT2[w][y][r]-np.mean(KOthenWT2[w][y]))))
            r+=1
        wtsd = np.sqrt((sum(sdArray))/(len(KOthenWT2[w][y])-1))
        
        SdWT.append(wtsd)
        if wtmean != 0:
            wtcv = (wtsd/wtmean)*100
        else:
            wtcv =0
        CvWT.append(wtcv)
        y+=1
    x=0
    z=1
    w+=1
#formatting output file
ry=0
formattedFinal=[]
labelAr=[]
labelAr.append('Peptide')
labelAr.append('Transition')
labelAr.append('Rank')
labelAr.append('Modified Sequence')
labelAr.append('Precursor M/Z')
labelAr.append('Protein Name')
labelAr.append('Protein Description')
labelAr.append('Precursor Charge')
labelAr.append('Retention Time')
u=0
while u<len(condAr):
    labelAr.append('Mean '+condAr[u])
    labelAr.append('Stan Dev '+condAr[u])
    labelAr.append('CV '+condAr[u])
    u+=1
i=0
j=0
while ry<len(pepAr):
    formatAr=[]
    formatAr.append(pepAr[ry][1])
    formatAr.append(ionAr[ry][1])
    formatAr.append(rankAr[ry][1])
    formatAr.append(mseqAr[ry][1])
    formatAr.append(precMzAr[ry][1])
    formatAr.append(proNameAr[ry][1])
    formatAr.append(prodesAr[ry][1])
    formatAr.append(precChAr[ry][1])
    formatAr.append(rtAr[ry][1])
    x=0
    while x < conditions:
        formatAr.append(MnWT[i])
        formatAr.append(SdWT[i])
        formatAr.append(CvWT[i])
        i+=1  
        x+=1
    formatAr.append('')
    ry+=1
    formattedFinal.append(formatAr)

#outputting csv statistics 
fn = os.path.splitext(fn)[0] + 'OUTPUT.csv'
with open(fn, 'w') as myfile:
    outputFile = csv.writer(myfile)
    i=0
    outputFile.writerow(labelAr)
    while i<len(formattedFinal):
        outputFile.writerows([formattedFinal[i]])
        i+=1
myfile.close()
