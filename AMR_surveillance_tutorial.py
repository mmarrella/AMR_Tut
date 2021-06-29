#!/usr/bin/env python
# coding: utf-8

# ## Tutorial Outline
# 
#     - Downloading data from NCBI FTP server
#     - 
#   
# 

# In[1]:


import ftplib
import os
from ftplib import error_perm
import re
from numpy import split
import pandas as pd
import numpy as np
import XlsxWriter


# In[2]:


def int_check(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    # Human ssorting is required to remove error in general python sorting which sorts only on the basis iof first digit which makes 854 > 2700
    # http://nedbatchelder.com/blog/200712/human_sorting.html
    return [ int_check(c) for c in re.split(r'(\d+)', text) ]



# In[5]:


ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
ftp.login()
ftp.cwd("/pathogen/Results/Escherichia_coli_Shigella")
# list files with ftplib
file_list = ftp.nlst() #returns the names of the elements present in a directory
filelist=[]
print(file_list)
for f in file_list: 
    #print (f)
    try:
        if 'pdg000000004' in f.lower():
            ftp.cwd("/pathogen/Results/Escherichia_coli_Shigella/" +f+"/AMR")
            file_list = ftp.nlst()
            #print(file_list)
            
            for f1 in file_list:
                if "pdg000000004" in f1.lower() and any(f1.endswith(ext) for ext in ['tsv']):
                    #print(f1)
                    filelist.append(f1)
                    #print("if loop ends here")
                    
            
    except ftplib.error_perm as e:
        print("AMR does not exist")
  
ftp.quit()




# In[6]:



filelist.sort(key=natural_keys)
f1_d=filelist[-1]
print(f1_d)



# In[16]:


#%%
ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
ftp.login()
ftp.cwd("/pathogen/Results/Escherichia_coli_Shigella")
# list files with ftplib
file_list = ftp.nlst()
#print(file_list)
for f in file_list: 
    #print (f)
    try:
        if 'pdg000000004' in f.lower():
            ftp.cwd("/pathogen/Results/Escherichia_coli_Shigella/" +f+"/AMR")
            file_list = ftp.nlst()
            #print(file_list)
            
            for f1 in file_list:
                if f1_d == f1 :
                    print(f1)
                    ftp.retrbinary("RETR {}".format(f1), open(f1, "wb").write)
                    print("downloaded {}".format(f1))       
              
    except ftplib.error_perm as e:
        continue
ftp.quit()
    


# In[59]:


import pandas as pd
df_Ecoli= pd.read_csv('/Users/jhasneha/Downloads/isolates.tsv', sep='\t', dtype='unicode')


# In[61]:



# dataframe.ndim
shape = df_Ecoli.shape
print("Dimension of the data: ",shape ) 


# In[62]:


df_Ecoli.columns
#df.head(10)


#%%


df_Ecoli["Scientific name"].value_counts(sort=True)
df_EC = df_Ecoli[df_Ecoli["Scientific name"] =="Escherichia coli" ]
df_EC["Scientific name"].value_counts(sort=True)


# In[66]:


df_EC["Location"].value_counts(sort=True)
df_EC_US = df_EC.loc[df_EC["Location"].str.contains('USA') == True]
df_EC_US["Location"].value_counts(sort=True)


# In[67]:


df_EC_US["Isolate"].value_counts(sort=True)


# In[68]:


df_EC_US["bioproject_center"].value_counts(sort=True)


# In[69]:


df_EC_US["isolation_source"].value_counts(sort=True)


# In[70]:


df_EC_US["host"].value_counts(sort=True)


# In[71]:


df_EC_US["source_type"].value_counts(sort=True)


# In[72]:


df_EC_US["epi_type"].value_counts(sort=True)

#%%
df_EC_US["AMR genotypes"].value_counts(sort=True)



# In[82]:


df = df_EC_US[["AMR genotypes"]].copy()
dfhead=df.head(100)
#print ("the pandas dataframe ",df
#%%
def create_df(df):
    df.columns=["AMR_genotypes"]
    df=df.replace(r'^\s*$','VALUE=ND', regex=True)
    df= df.fillna(0)
    df=df.replace(0, 'VALUE=ND', regex=True)
    df['AMRgenotypes']=df['AMR_genotypes']
    df=df.set_index('AMRgenotypes')
    split_df= df['AMR_genotypes'].str.split(',', expand= True)
    print (split_df.head())
    #now we make the column list for the table using unique values
    split_df_list=[]
    split_df_list = split_df.values.tolist()
    flat_list_sd = []
    for sublist in split_df_list:
        for item in sublist:
            flat_list_sd.append(item)
    split_flsd=[]
    for element in flat_list_sd:
        if element != None:
            split_flsd.append(element.split('=', 1)[0]) 
    header_set=set(split_flsd)
    header_list_AST=list(header_set)
    df_AST=pd.DataFrame(header_list_AST)
    df_AST.columns=["AMR"]
    df_AST=df_AST.T
    df_AST=df_AST.reset_index()
    df_AST=df_AST.T
    df_AST=df_AST.reset_index()
    df_AST=df_AST.drop(['index'], axis=1)
    df_AST.columns=["AMR"]
    print("the columns are,",df_AST)
    x=1
    for tup in split_df.itertuples():
        try:
            #print (tup)
            tuplist=list(tup)
            tup_list_df=pd.DataFrame(tuplist)
            colname=tup_list_df.iloc[0]
            head= pd.DataFrame([colname])
            head=head.reset_index()
            head= head.T
            head=head.reset_index()
            head=head.drop(columns=['index'])
            #print(head)
            head.at[0,0]="AMR"
            head.columns=["header"]
            head_list=head['header'].tolist()
            tup_list_df=tup_list_df.drop([0])
            tup_list_df.columns=["AMR"]
            #print(tup_list_df)
            tup_list_df= tup_list_df.AMR.str.split('=', expand= True)
            tup_list_df.columns=head_list
            
            tup_list_df=tup_list_df.T
            tup_list_df=tup_list_df.reset_index()
            tup_list_df=tup_list_df.T
            tup_list_df.columns=["AMR","AMR_g"]
            #print(tup_list_df)
            #tup_list_df_list= list(tup_list_df)
            tup_list_df_list=list(tup_list_df.itertuples(index=False, name=None))
            tup_list_df_list_dict=dict(tup_list_df_list)
            #tup_list_df_list_dict=dict(tup_list_df_list)
            df_AST[x]=df_AST['AMR'].map(tup_list_df_list_dict)
            x=x+1
            print (x)

        except:
            print("No result updated")
            continue
    
    df_AST=df_AST.T
    return (df_AST)

#%%
df_AMRgene=create_df(df)
print("AMR genotype dataframe created")

#%%

df_AMRgene.columns = df_AMRgene.iloc[0] #grab the first row for the header
df_AMRgene=df_AMRgene.rename(columns={'AMR':'AMR_genotype'})
df_AMRgene=df_AMRgene.reset_index()# the change of index in the previous command causes error in the concatenation of the other dataframe
df_AMRgene= df_AMRgene.drop(columns=['index'])# so reset is used to create a new index and the old index is deleted 
df_AMRgene=df_AMRgene.drop([0])# drops the first row which is a header row from dataframe but this also changes the index of the dataframe
df_AMRgene= df_AMRgene.drop(columns=['VALUE'])
df_AMRgene=df_AMRgene.reset_index()
df_AMRgene= df_AMRgene.drop(columns=['index'])
df_AMRgene=df_AMRgene.reset_index()
df_AMRgene=df_AMRgene.rename(columns={'index':'AMR_id'})

#%%
df_AMRgene.to_csv("df_AMEGENE.tsv", sep='\t')
df_AMRgene.to_excel('df_AMRGENE.xlsx', sheet_name='Sheet1', engine='xlsxwriter')



# %%
df_test = df_EC_US[["AST phenotypes"]].copy()

#dfhead=df.head(20)
# %%
df_astresult=create_df(df_test)


# %%
df_astresult.columns = df_astresult.iloc[0] #grab the first row for the header
df_astresult=df_astresult.rename(columns={'AMR':'AST_phenotype'})
df_astresult=df_astresult.reset_index()# the change of index in the previous command causes error in the concatenation of the other dataframe
df_astresult= df_astresult.drop(columns=['index'])# so reset is used to create a new index and the old index is deleted 
df_astresult=df_astresult.drop([0])# drops the first row which is a header row from dataframe but this also changes the index of the dataframe
df_astresult= df_astresult.drop(columns=['VALUE'])
df_astresult=df_astresult.reset_index()
df_astresult= df_astresult.drop(columns=['index'])
df_astresult=df_astresult.reset_index()
df_astresult=df_astresult.rename(columns={'index':'AST_id'})

#%%
df_astresult.to_csv("df_astresult.tsv", sep='\t')



#%%
df_EC_US=df_EC_US.reset_index()
df_EC_US=df_EC_US.rename(columns={'index':'ID'})
# %%
df_main=pd.DataFrame()
df_main=df_EC_US.merge(df_astresult,how='left', left_on='ID', right_on='AST_id')
df_main=df_main.merge(df_AMRgene,how='left', left_on='ID', right_on='AMR_id')

# %%
df_Ecoli_SH=pd.DataFrame()
df_Ecoli_SH=df_EC_US[["Isolation Source","Host","Isolation type"]].copy()
df_Ecoli_SH["sourcehostcomb"] = df_EC_US["Isolation Source"].fillna('IMV').astype(str) +'_'+ df_EC_US["Host"].fillna('HMV').astype(str) +'_'+df_EC_US["Isolation type"].fillna('TMV').astype(str)
df_Ecoli_SH=df_Ecoli_SH.rename(columns={'Isolation Source':'Isolation_Source','Isolation type':'Isolation_type'})
# %%
Ecoli_IB_Merged_df=pd.crosstab( df_Ecoli_SH.sourcehostcomb,df_Ecoli_SH.Isolation_type,  margins=True, margins_name="total")

# %%
import re
Ecoli_IB_Merged_df=Ecoli_IB_Merged_df.reset_index()
Ecoli_IB_Merged_df["epitype"]=Ecoli_IB_Merged_df["sourcehostcomb"]
Ecoli_IB_Merged_df['epitype'] = Ecoli_IB_Merged_df['epitype'].map(lambda x: re.sub(r'W+', '-', x)) #removes all characters which are not words 
Ecoli_IB_Merged_df['epitype'] = Ecoli_IB_Merged_df['epitype'].map(lambda x: re.sub(" ", "-", x))
Ecoli_IB_Merged_df['epitype'] = Ecoli_IB_Merged_df['epitype'].str.lower()#converts the colmn into all lower cases
Ecoli_IB_Merged_df["epi_hostname"]=Ecoli_IB_Merged_df["epitype"]
Ecoli_IB_Merged_df.loc[Ecoli_IB_Merged_df.epi_hostname.str.contains('clinical|sapiens',na=False),'epi_hostname'] = 'Homosapiens'
Ecoli_IB_Merged_df.loc[Ecoli_IB_Merged_df.epi_hostname.str.contains('stockpellet|livestockfeed|animalfeed|animalfeed|milksus|milkporcine|petfood|peaproteinpellet|animalfeed|grass|feedgrass|feedhmv|catfood|dogfood',na=False),'epi_hostname'] = 'Animalfood'
Ecoli_IB_Merged_df.loc[Ecoli_IB_Merged_df.epi_hostname.str.contains('cheese|chocolate|icecream|dairyproducts|milk|roquefort',na=False),'epi_hostname'] = 'Dairy'
Ecoli_IB_Merged_df.loc[Ecoli_IB_Merged_df.epi_hostname.str.contains('veggie|piecrust|lettuce|spinach|coriander|sprouts|pepper|flour|alfalfa|leafygreens|basil|cherrytomato|cucumber|parsley|pizzadough|kale|grain|celery|almond|berries|berry|applecider|apple|cider|chard|cilantro|clover|soynut|tomato|carrot|forage|mint|springmix|vegetable|romaine|melon|microgreen|nuts|organic|oregano|basil|freshsage|rosemary|bean|foodhmv|foodsamplehmv|foodsamplefood|cepa|kolarabi|peaprotein|cabbage|leaf|thyme|springmix|creamysoynutbutter|soynut|cantaloupehmv|cantalopehmv',na=False),'epi_hostname'] = 'plant_based_food'
Ecoli_IB_Merged_df.loc[Ecoli_IB_Merged_df.epi_hostname.str.contains('porkhmv|wholechicken|chickenskin|skinsamplegallus|foodturkey|egg|vealtrim|beefstreak|carcass|carcassswab|beefcarctrim|beefcutlet|beefpatties|patties|turkeypatties|rawveal|rawbeef|beeftrim|groundbeef|groundpork|groundmeat|groundturkey|meatsticks|meatsauce|tenderloin|hamburger|burger|meat|meatloaf|commercial|salami|shrimp|piecrust|venison|boneless|breast|steak|siluriformes|porkchop|product|ground|comminuted|food(pork)|foodpork|packaged|sliced|yolk|chickenwing|chickenthigh|wing|chickenleg|mixedparts|beefcarcass|venison',na=False),'epi_hostname'] = 'Meat/Meatproduct'
Ecoli_IB_Merged_df.loc[Ecoli_IB_Merged_df.epi_hostname.str.contains('bullmanure|ovine|bostaurus|beef|cattle|cow|bos|bovine|calf|veal|buffalo|goat|udder|mastitis|sheep|caprahircus|ovisaries',na=False),'epi_hostname'] = 'Ruminant'
Ecoli_IB_Merged_df.loc[Ecoli_IB_Merged_df.epi_hostname.str.contains('pig|porcine|swine|sus|feralswine|feralpig',na=False),'epi_hostname'] = 'Swine'
Ecoli_IB_Merged_df.loc[Ecoli_IB_Merged_df.epi_hostname.str.contains('turkey|chicken|goose|duck|gallus|meleagris',na=False),'epi_hostname'] = 'Poultry'
Ecoli_IB_Merged_df.loc[Ecoli_IB_Merged_df.epi_hostname.str.contains('canine|dog|familiaris|familaris',na=False),'epi_hostname'] = 'Dog'
Ecoli_IB_Merged_df.loc[Ecoli_IB_Merged_df.epi_hostname.str.contains('horse|equus|columba|equine|equis',na=False),'epi_hostname'] = 'Horse'
Ecoli_IB_Merged_df.loc[Ecoli_IB_Merged_df.epi_hostname.str.contains('fecesfeline|lungcat|catfeline|imvcat|cathmv|feliscatus|urinefeline',na=False),'epi_hostname'] = 'Cat'
Ecoli_IB_Merged_df.loc[Ecoli_IB_Merged_df.epi_hostname.str.contains('environment|environmental',na=False),'epi_hostname'] = 'Environment'

#%%
Ecoli_IB_Merged_df=Ecoli_IB_Merged_df.rename(columns={'epi_hostname':'OneHealth_host'})
#%%
Ecoli_IB_Merged_df.to_csv("Ecoli_IB_Merged_df.tsv", sep='\t')

# %%
df_main["sourcehostcomb"] = df_main["Isolation Source"].fillna('IMV').astype(str) +'_'+ df_main["Host"].fillna('HMV').astype(str) +'_'+df_main["Isolation type"].fillna('TMV').astype(str)
Oh_dict = dict(zip(Ecoli_IB_Merged_df.sourcehostcomb, Ecoli_IB_Merged_df.OneHealth_host))
df_main['OneHealth_host']=df_main['sourcehostcomb'].map(Oh_dict)

# %%
df_main["OneHealth_host"].value_counts()

# %%

df_main["Collection Date"].value_counts()

#%%
df_CD=df_main[["Collection Date"]].copy()
def cd(df_CD):
    df_CD = df_CD.fillna('2030') 
    df_CD = df_CD.replace('1900/1949','1949') 
    df_CD = df_CD.replace('1949/1900','1949')
    df_CD = df_CD.replace('1900/1960','1960')
    df_CD = df_CD.replace('2017/2018','2018')
    df_CD = df_CD.replace('2015-01/2016-07','2015')
    df_CD = df_CD.replace('2018-11/2019-03','2019')
    df_ser=pd.Series(df_CD.values.tolist(), index=df_CD.index)
    df_CD['cd_mod']=df_ser.apply(lambda x: pd.to_datetime(x).strftime('%m-%d-%Y')[0])
    df_CD['cd_year']=df_ser.apply(lambda x: pd.to_datetime(x).strftime('%Y')[0])
    return (df_CD)

df_date=cd(df_CD)
df_date=df_date.reset_index()
df_date=df_date.rename(columns={'index':'cd_id'})
df_main=df_main.merge(df_date,how='left', left_on='ID', right_on='cd_id')
df_main = df_main["cd_year"].replace('2030',np.nan) 
# %%
df_main["Location"].value_counts()
#%%
gloc_df=df_main[["Location"]].copy()


# %%
df_stcode=pd.read_csv("statelatlong.csv")
def locdict(gloc_df, df_stcode):
    loc_dict=df_stcode.set_index('statecode')['state'].to_dict()
    stcode_list=df_stcode['statecode'].tolist()
    st_list=df_stcode['state'].tolist()
    gloc_df["state"]=gloc_df["Location"]
    gloc_df["state"]=gloc_df.state.str.split(':').str[1]
    gloc_df["state"] = gloc_df.state.str.lstrip() 
    #extracting the state code by using pattern 'pat' and assigns na where not available
    pat = r'({})'.format('|'.join(stcode_list))
    gloc_df['state'] = gloc_df['state'].str.extract(pat, expand=False).fillna(gloc_df['state'])
    gloc_df['state'] =gloc_df['state'].map(loc_dict).fillna(gloc_df['state'])
    #the same as statecode but the pattern list is state name
    pat = r'({})'.format('|'.join(st_list)) 
    gloc_df['state'] = gloc_df['state'].str.extract(pat, expand=False).fillna(gloc_df['state'])
    #replacing the center names
    gloc_df=gloc_df.replace(to_replace ="Fort Sam Houston", value ="Texas")
    gloc_df=gloc_df.replace(to_replace ="Houston", value ="Texas")
    gloc_df=gloc_df.replace(to_replace ="Baltimore", value ="Maryland")
    gloc_df=gloc_df.replace(to_replace ="St. Louis Clyde watershed of Lake Superior", value ="St Louis")
    gloc_df=gloc_df.replace(to_replace ="UC Davis Medical Center", value ="California")
    gloc_df=gloc_df.replace(to_replace ="UC Davis Medical Center, Davis, Ca", value ="California")
    gloc_df=gloc_df.replace(to_replace ="University of Miami Department of Pathology", value ="Florida")
    gloc_df=gloc_df.replace(to_replace ="Burlington", value ="Vermont")
    gloc_df=gloc_df.replace(to_replace ="Lahey Hospital & Medical Center", value ="Massachusetts")
    gloc_df=gloc_df.replace(to_replace ="Ronald Reagan UCLA Medical Center", value ="California")
    gloc_df=gloc_df.replace(to_replace ="Santa Barbara, UCSB animal house", value ="New Jersey")
    gloc_df=gloc_df.replace(to_replace ="Robert Wood Johnson", value ="California")
    gloc_df=gloc_df.replace(to_replace ="ALACHUA", value ="Florida")
    gloc_df["state"]=gloc_df["state"].replace ("USA",np.nan)

    return ( gloc_df)
df_state=locdict(gloc_df,df_stcode)
# gdf=df_state.copy()
# gdf=gdf.drop_duplicates(subset ="Location",
#                         keep = False, inplace = True)
#     #write it to a file to be used as dictionary for further pathogens 
# gdf.to_csv("Gloc_updated.csv", index=False)

#%%
st_dict = dict(zip(df_state.Location, df_state.state))
df_main['Collection_state']=df_main['Location'].map(st_dict)


# %%
df_main['Collection_state'].value_counts()


# %%
df_main.to_csv("df_main_complete.tsv", sep="\t", index=False )

""".............
Analysis of the data
........................."""

# %%
df_mainC=pd.read_csv("df_main_complete.tsv", sep="\t",dtype='unicode')
#%%
df_astresult=pd.read_csv("df_astresult.tsv", sep="\t")
list_ast=list(df_astresult.columns)
list_ast.remove('AST_phenotype')

# %%
df_AMRgene=pd.read_excel("df_AMRgene.xlsx")
df_AMRgene = df_AMRgene.iloc[: , 1:] #drops the first index column
list_amrgene=list(df_AMRgene.columns)
list_amrgene.remove('AMR_genotype')



# %%
def gene_frequency(lst_,df_kleb):
    try:
        appended_data=[]
        for k in lst_:
            try:
                xdf_crosstab=pd.crosstab(df_kleb["Isolate"],df_kleb[k])
                xdf_crosstab["AMR_gene"] = k
                print('xdf',xdf_crosstab)
                print('') # for spacing
                if len(xdf_crosstab.index) >1:
                    #print("xdf_crosstab.index",xdf_crosstab.index)
                    appended_data.append(xdf_crosstab)
                else:
                    print("only one tested")
                    
            except:
                print("not found")
                continue
        appended_data= pd.concat(appended_data, sort= False) 
        appended_data=appended_data.drop("missing")
        appended_data=appended_data.drop("")
    except:
        print("completed??",k)

    return(appended_data)
#%%
df_gf=gene_frequency(list_amrgene,df_mainC)
df_gf = df_gf.replace( np.nan, 0) 
df_gf["Gene_frequency"]= df_gf.sum(axis=1)
df_gf= df_gf.drop(["MISTRANSLATION","PARTIAL","POINT","PARTIAL_END_OF_CONTIG","COMPLETE","HMM"], axis=1)
    
#%%
df_gc=df_gf.groupby(['AMR_gene']).sum()
df_gc=df_gc.reset_index()


# creating the ECDF plot for selecting Important genes
#%%
from matplotlib import pyplot
import matplotlib.pyplot as plt
from statsmodels.distributions.empirical_distribution import ECDF

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 600
sample1 = np.array(df_gc["Gene_frequency"])
#sample2 = normal(loc=40, scale=5, size=700)
#sample = hstack((sample1, sample2))
# fit a cdf
ecdf = ECDF(sample1)

# get cumulative probability for values
print('P(x<100): %.3f' % ecdf(100))
print('P(x<500): %.3f' % ecdf(500))
print('P(x<1200): %.3f' % ecdf(1200))
# plot the cdf
pyplot.plot(ecdf.x, ecdf.y)
plt.title("AMR Gene frequency ECDF plot for $\it{E.coli}$ AMR Relational data model",fontsize=10)
plt.ylabel("Probability of occurence of AMR gene in the data model",fontsize=10)
plt.xlabel("gene frequency",fontsize=8)
plt.axhline(y=0.9, color='m', linewidth= 0.8)
#plt.savefig("/Users/jha/Documents/Spring2021/SPR_indot/graphs/Recdf12.pdf", dpi=900)
pyplot.show()
plt.title("AMR Gene frequency plot for $\it{E.coli}$ AMR Relational data model",fontsize=8)
plt.ylabel("Number of AMR genes in the data model",fontsize=8)
plt.xlabel("gene frequency", fontsize=10)
pyplot.hist(sample1)
#plt.savefig("/Users/jha/Documents/Spring2021/SPR_indot/graphs/Rhist12.pdf", dpi=900)

pyplot.show()
# %%
ecdf_lim=1200
df_gc = df_gc[df_gc["Gene_frequency"] >=ecdf_lim] 
list_impgenes=list(df_gc["AMR_gene"])
df_gc=df_gc.set_index("AMR_gene")
# %%

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 600
ax=df_gc.plot.bar(color = "blue") #the dataframe object plots according to its index
plt.legend(loc='upper right')
plt.ylabel('Gene_Frequency', fontsize=6)
plt.xlabel("AMR_gene",fontsize=6)
plt.xticks(fontsize=6)
plt.title("Gene frequency plot of the important genes", fontsize=8)

# %%
df_tf=gene_frequency(list_ast,df_mainC)
df_tf= df_tf.drop(["ND","SDD"], axis=1)
df_tf = df_tf.replace(np.nan, 0)  
df_tf["Non-Susceptible"]= df_tf["I"]+df_tf["R"]
df_tf= df_tf.drop(["I","R"], axis=1)

# %%
df_tc=df_tf.groupby(['AMR_gene']).sum()
df_tc=df_tc.reset_index() #flattens the column name
#%%
df_tc=df_tc.rename(columns={'AMR_gene':'Antimicrobial'})
df_tc=df_tc.rename(columns={'S':'Susceptible'})
df_tc=df_tc.set_index("Antimicrobial")

# %%
def plot_ast(df,string, string2):
    import matplotlib as mpl
    mpl.rcParams['figure.dpi'] = 600
    import matplotlib
    from matplotlib.cm import viridis
    ax=df.plot.bar(stacked=True)#the dataframe object plots according to its index
    plt.legend(bbox_to_anchor=(1,1), loc='upper left')
    ax.set_ylabel('Number of isolates', fontsize=8)
    ax.set_xlabel(string2,fontsize=8)
    plt.xticks(fontsize=6, rotation=80)
    plt.title(string, fontsize=8)
    return (plt.show())

# %%

plot_ast(df_tc,"AST test result for list of Antimicrobials")
# %%
plot_ast(df_gc,"Gene frequency plot of the important genes")
#%%
# can use the funtion to plot the collection date, location and onehealth_host
# df_sh=df_mainC["OneHealth_host"].value_counts()
# plot_ast(df_sh,"Number of Isolates tested from each sample host","OneHealth Host")
# df_cd=df_mainC["cd_year"].value_counts()
# plot_ast(df_cd,"Number of Isolates tested from each year","Collection Year")


"""================================================


Creating the table with all values for the analysis for important genes

==================================================="""


# %%
df_imp=gene_frequency(list_impgenes,df_mainC)
df_imp=df_imp.reset_index()
df_imp = df_imp.replace( np.nan, 0) 
df_imp["Gene_frequency"]= df_imp.sum(axis=1)
df_imp= df_imp.drop(["MISTRANSLATION","PARTIAL","POINT","PARTIAL_END_OF_CONTIG","COMPLETE","Gene_frequency"], axis=1)

# %%
def create_all(df_kleb,appended_data):
    isoyear_dict = dict(zip(df_kleb.Isolate, df_kleb.cd_year))
    appended_data['CD_year']=appended_data['Isolate'].map(isoyear_dict)
    isostate_dict = dict(zip(df_kleb.Isolate, df_kleb.Collection_state))
    appended_data['Collection_state']=appended_data['Isolate'].map(isostate_dict)
    isosamp_dict = dict(zip(df_kleb.Isolate, df_kleb.OneHealth_host))
    appended_data['OH_host']=appended_data['Isolate'].map(isosamp_dict)
    return(appended_data)
df_imp_all=create_all(df_mainC,df_imp)

#%%
df_imp_all.head()


# %%
dfgsh = pd.pivot_table(df_imp_all,index=["OH_host"],columns=["AMR_gene"],values=["Isolate"],aggfunc=pd.Series.nunique, fill_value= 'None')
dfgsh=dfgsh.reset_index()
dfgsh.columns = dfgsh.columns.droplevel()
dfgsh=dfgsh.rename(columns={dfgsh.columns[0]:'OH_host'})
dfgsh=dfgsh.drop(columns=["acrF","blaEC","mdtM"])
dfgsh=dfgsh.T
dfgsh.columns=dfgsh.iloc[0]
dfgsh=dfgsh.reset_index()
dfgsh=dfgsh.drop([0])
dfgsh=dfgsh.set_index("AMR_gene")
dfgsh.head()


# %%
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 600
import matplotlib
from matplotlib.cm import tab10
n=len(list(dfgsh.columns))
def get_n_colors(n):
    return[ tab10(float(i)/n) for i in range(n) ]
colors = get_n_colors(n)
#viridis = cm.get_cmap('viridis', 14)
ax=dfgsh.plot.bar(stacked=True,color = colors)
plt.legend(bbox_to_anchor=(1,1), loc='upper left')
ax.set_ylabel('AMR genes', fontsize=8)
ax.set_xlabel('Number of Isolates',fontsize=8)


