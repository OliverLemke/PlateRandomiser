#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 11:33:53 2020

@author: Oliver
"""

import numpy as np
import pandas as pd
from sklearn.utils import shuffle
import string
import Color_Mix
import matplotlib.pyplot as plt
from matplotlib import colors

### Try to include multiple samples for one patient into the script
### Design Icon
### Automatically adjust Size for bigger window?

### Input ###
Plate_layout = "96" # switch between two plates 96 an 385
num_Blanks = 1 # default
# References
Reference_1 = "Plasma"
Reference_2 = "Serum"
### Either
num_References = 12
percentage_Ref_1 = 0.5 # Default None
### Or
num_Ref_1 = 6 # Default 6 or num_references
num_Ref_2 = 6 # Default 6
# minimal number of Plates
num_plates = 2
## Input file
in_file = "./Data.xlsx"
## Ouput file
out_file = "./Layout.txt"
## Columns for Properties (to be distributed equally) max. 3-4
imp_columns = ["Geschlecht","PT_Tumorgruppe","PT_UICC-Stadium","HPV"]
## Reference column (Label)
ref_columns = ["Interne Probennr","MN_Patient"]
## Sheet names (to be read in)
sheet_names = ["Tabelle1"]
## Seed for randomization steps
seed = 40

# Select plate Layout
def get_Plate_layout(Plate_layout, num_Blanks, num_References=None, num_Ref_1=6, num_Ref_2=6):
    """
    Get the plate layout for 96- or 384-well plate.

    Parameters
    ----------
    Plate_layout : {"96","384"}
        Number of wells per plate.
    num_Blanks : int
        Number of Blanks to be included per plate.
    num_References : int, optional
        Minimal number of references (if provided, for fixed numbers define num_Ref_1 and num_Ref_2). The default is None.
    num_Ref_1 : int, optional
        Number of Reference 1 (needed if num_References is not defined). The default is 6.
    num_Ref_2 : int, optional
        Number of Reference 2 (needed if num_References is not defined). The default is 6.

    Returns
    -------
    num_columns : int
        Number of columns per plate.
    num_rows : int
        Number of rows per plate.
    num_wells : int
        Number of wells per plate.
    num_wells_to_fill : int
        Maximal number of wells to be filled per plate.
    """

    if not num_References:
        num_References = num_Ref_1 + num_Ref_2
    
    if Plate_layout in ["96","384"]:
        if Plate_layout == "96":
            num_columns = 12
            num_rows = 8
        elif Plate_layout == "384":
            num_columns = 24
            num_rows = 16
        num_wells = num_columns*num_rows
        num_wells_to_fill = num_wells - num_Blanks - num_References
    else:
        raise ValueError("Plate layout unknown")
    return num_columns, num_rows, num_wells, num_wells_to_fill

# Set seed
def set_seed(seed=-1):
    """
    Set seed for randomization.

    Parameters
    ----------
    seed : int, optional
        Seed for randomization. The default is -1.
    """

    np.random.seed(seed)
    return

# Read data: add tsv, csv and xlsx
def read_data_excel(in_file, ref_columns, imp_columns, sheet_names=["Tabelle 1"]):
    """
    Read in Excel data.

    Parameters
    ----------
    in_file : str
        Name of the input file.
    ref_columns : list of str
        Identifier columns for the sample.
    imp_columns : list of str
        Features used for Fingerprint determination.
    sheet_names : list of str, optional
        Sheets to be used. The default is ["Tabelle 1"].

    Returns
    -------
    data_dict : dict
        Dictionary containing the data with ref_columns[0] as key.
    Fingerprints : array
        Fingerprints for all samples (row-wise).

    """
    data = pd.read_excel(in_file, sheet_name=sheet_names)
    data_merged = pd.concat(data)
    ## Save to dictionary ID: Info and extract important columns
    data_dict = {}
    Fingerprints = []
    for index,row in data_merged.iterrows():
        data_dict.update({row[ref_columns[0]]:row})
        Fingerprints.append(row[imp_columns])
    return data_dict, Fingerprints

# Get Fingerprint
def get_Fingerprint(data_dict, Fingerprints, imp_columns):
    """
    Group data according to their fingerprints and assign samples to common fingerprint.

    Parameters
    ----------
    data_dict : dict
        Dictionary containing the data with ref_columns[0] as key.
    Fingerprints : array
        Fingerprints for all samples (row-wise).
    imp_columns : list of str
        Features used for Fingerprint determination.

    Returns
    -------
    Fingerprint_IDs : array
        Unique Fingerprints.
    Fingerprints_list : list of str
        Samples assigned to each Fingerprint.
    """
    
    Fingerprints = np.asarray(Fingerprints,dtype="U25")
    Fingerprint_IDs = np.unique(Fingerprints,axis=0)
    Fingerprint_dict = {}
    for index,fingerprint in enumerate(Fingerprint_IDs):
        Fingerprint_dict.update({index:fingerprint})
    
    Fingerprints_list = []
    for fingerprint_index in Fingerprint_dict:
        a = []
        for ind,key in enumerate(data_dict):
            if all(data_dict[key][imp_columns].values.astype("U25")==Fingerprint_dict[fingerprint_index]):
                a.append(key)
        
        Fingerprints_list.append(shuffle(np.asarray(a)))
    return Fingerprint_IDs, Fingerprints_list

# Get Fingerprint statistics

def get_Fingerprint_statistics(Fingerprints_list, num_wells_to_fill, num_plates):
    """
    Get statistics for every fingerprint and determine how many samples are added per plate.

    Parameters
    ----------
    Fingerprints_list : list of str
        Samples assigned to each Fingerprint.
    num_wells_to_fill : int
        Maximal number of wells to be filled per plate.
    num_plates : int
        Minimal number of plates.

    Returns
    -------
    num_Fingerprint : array
        Number of samples per fingerprint
    num_plates : int
        Final number of plates.
    num_Sample_per_plate : array
        Number of samples per Fingerprint per plate.
    num_Samples_overlap : array
        Number of samples not assigned to any plate.
    """
    
    num_Fingerprint = np.asarray([len(item) for item in Fingerprints_list])
    num_plates_sug = int(np.ceil(np.sum(num_Fingerprint)/num_wells_to_fill))
    if num_plates<num_plates_sug:
        num_plates = num_plates_sug
        print("Number of Plates changed to %d" %num_plates)
    
    ## Determine how many per plate (floor)
    num_Sample_per_plate = np.int_(np.floor(num_Fingerprint/num_plates))
    
    ## Determine how many are too large (mod)
    num_Samples_overlap = np.int_(np.mod(num_Fingerprint,num_plates))
    return num_Fingerprint, num_plates, num_Sample_per_plate, num_Samples_overlap

def distribute_Samples(Fingerprints_list, data_dict, num_Sample_per_plate, num_plates):
    """
    Distribute Samples across the plates

    Parameters
    ----------
    Fingerprints_list : list of str
        Samples assigned to each Fingerprint.
    data_dict : dict
        Dictionary containing the data with ref_columns[0] as key.
    num_Sample_per_plate : array
        Number of samples per Fingerprint per plate.
    num_plates : int
        Final number of plates.

    Returns
    -------
    Plates_overlap : list of arrays
        Samples distributed over num_plates plates.
    """
    
    Plates = []
    for i in range(num_plates):
        plate = np.asarray([])
        for ind,fingerprint in enumerate(Fingerprints_list):
            ## Cut at pre-determined cut
            plate = np.concatenate((plate,fingerprint[i*num_Sample_per_plate[ind]:(i+1)*num_Sample_per_plate[ind]]))
        Plates.append(plate)
        
    ## Distribute remaining data points
    Overlap = np.asarray([])
    
    for i in range(len(Fingerprints_list)):
        Overlap = np.concatenate((Overlap,Fingerprints_list[i][num_Sample_per_plate[i]*num_plates:]))
    Overlap = shuffle(Overlap)
    Plates_overlap = []
    
    for ind,plate in enumerate(Plates):
        Plates_overlap.append(np.concatenate((plate,Overlap[ind::num_plates])))
        
    ## Check whether all samples are distributed (Compare Input and Output)
    Check = np.asarray([item for plate in Plates_overlap for item in plate])
    print("Unique IDs provided: %d" %len(data_dict))
    print("Unique IDs detected on plates: %d" %len(np.unique(Check)))
    return Plates_overlap

def distribute_References(Plates_overlap, Reference_1, Reference_2, num_wells, num_Blanks, num_Ref_1=6, num_Ref_2=6, percentage_Ref_1=None):
    """
    Add references and blanks to the plate

    Parameters
    ----------
    Plates_overlap : list of arrays
        Samples distributed num_plates over plates.
    Reference_1 : str
        Label of Reference 1.
    Reference_2 : str
        Label of Reference 2.
    num_wells : int
        Number of wells per plate.
    num_Blanks : int
        Number of blanks to be included per plate.
    num_Ref_1 : int, optional
        Number of Reference 1 (needed if num_References is not defined). The default is 6.
    num_Ref_2 : int, optional
        Number of Reference 2 (needed if num_References is not defined). The default is 6.
    percentage_Ref_1 : float, optional
        Percentage of Reference 1 (needed if num_Ref_1 and num_Ref_2 are not defined). The default is None.

    Returns
    -------
    Plates_final : list of arrays
        Samples, References and Blanks distributed over num_plates plates.
    """

    Plates_final = []
    for ind,plate in enumerate(Plates_overlap):
        if percentage_Ref_1:
            num_Ref = int(num_wells-len(plate)-num_Blanks)
            num_Ref_1 = int(np.floor(percentage_Ref_1 * num_Ref))
            num_Ref_2 = int(np.floor((1-percentage_Ref_1) * num_Ref))
            if num_Ref_1+num_Ref_2 != num_Ref:
                p = np.random.uniform(0,1)
                if p<percentage_Ref_1:
                    num_Ref_1+=1
                else:
                    num_Ref_2+=1
        else:
            num_Blanks += int(num_wells-len(plate)-num_Blanks-num_Ref_1-num_Ref_2)
        Ref_1 = np.asarray([Reference_1]*num_Ref_1)
        Ref_2 = np.asarray([Reference_2]*num_Ref_2)
        Blanks = np.asarray(["Blank"]*num_Blanks)    
        Plates_final.append(shuffle(np.concatenate((plate,Ref_1,Ref_2,Blanks))))
    return Plates_final

def generate_Output(data_dict, Plates_final, ref_columns, imp_columns, num_columns, out_file="Out.txt"):
    """
    Write Output file (tab-separated)

    Parameters
    ----------
    data_dict : dict
        Dictionary containing the data with ref_columns[0] as key.
    Plates_final : list of arrays
        Samples, References and Blanks distributed over num_plates plates.
    ref_columns : list of str
        Identifier columns for the sample.
    imp_columns : list of str
        Features used for Fingerprint determination.
    num_columns : int
        Number of columns per plate.
    out_file : str, optional
        Name of the output file. The default is "Out.txt".
    """

    dict_alph, dict_num = get_Dictionaries()
    
    remaining_columns = [column for column in data_dict[list(data_dict.keys())[0]].index if (column not in ref_columns) & (column not in imp_columns)]
    
    print("Writing output")
    c_Ref_1 = 1
    c_Ref_2 = 1
    c_Blank = 1
    
    with open(out_file,"w") as file:
        for ref_column in ref_columns:
            file.write("%s\t" %ref_column)
        for column in imp_columns:
            file.write("%s\t" %column)
        file.write("Analysis.Plate\tAnalysis.Row\tAnalysis.Column\t")
        for ind,column in enumerate(remaining_columns):
            if ind != len(remaining_columns)-1:
                file.write("%s\t" %column)
            else:
                file.write("%s\n" %column)    
        for index,plate in enumerate(Plates_final):
            for index2,el in enumerate(plate):
                if el == Reference_1:
                    for i in range(len(ref_columns)):
                        file.write(Reference_1+"_{:03}".format(c_Ref_1)+"\t")
                    for i in range(len(imp_columns)):
                        file.write("QC_1\t")
                    c_Ref_1+=1
                elif el == Reference_2:
                    for i in range(len(ref_columns)):
                        file.write(Reference_2+"_{:03}".format(c_Ref_2)+"\t")
                    for i in range(len(imp_columns)):
                        file.write("QC_2\t")
                    c_Ref_2+=1
                elif el == "Blank":
                    for i in range(len(ref_columns)):
                        file.write("Blank_"+"{:03}".format(c_Blank)+"\t")
                    for i in range(len(imp_columns)):
                        file.write("Blank\t")
                    c_Blank+=1
                else:
                    for ref_column in ref_columns:
                        file.write("%s\t" %data_dict[el][ref_column])
                    for column in imp_columns:
                        file.write("%s\t" %data_dict[el][column])
                file.write(str(index+1)+"\t")
                file.write("%s\t" %dict_alph[np.floor(index2/num_columns)])
                file.write(str(np.mod(index2,num_columns)+1)+"\t")
                for ind,column in enumerate(remaining_columns):
                    if ind != len(remaining_columns)-1:
                        if (el != Reference_1) & (el != Reference_2) & (el != "Blank"):
                            file.write("%s\t" %data_dict[el][column])
                        else:
                            file.write("NA\t")
                    else:
                        if (el != Reference_1) & (el != Reference_2) & (el != "Blank"):
                            file.write("%s\n" %data_dict[el][column])
                        else:
                            file.write("NA\n")
    return

def get_Statistics(data_dict, Fingerprint_IDs, Plates_final, imp_columns):
    """
    Get statistics for the plate layout

    Parameters
    ----------
    data_dict : dict
        Dictionary containing the data with ref_columns[0] as key.
    Fingerprint_IDs : array
        Unique Fingerprints.
    Plates_final : list of arrays
        Samples, References and Blanks distributed over num_plates plates.
    imp_columns : list of str
        Features used for Fingerprint determination.

    Returns
    -------
    num_Fingerprint_per_plate : array
        Number of specific fingerprints on every plate.
    """

    num_Fingerprint_per_plate=np.zeros((len(Fingerprint_IDs),len(Plates_final)))
    for ind_p,plate in enumerate(Plates_final):
        for el in plate:
            for index,fingerprint in enumerate(Fingerprint_IDs):
                try:
                    if all(data_dict[el][imp_columns].values.astype("U25")==fingerprint):
                        num_Fingerprint_per_plate[index,ind_p]+=1
                except:
                    pass
    return num_Fingerprint_per_plate

def print_Statistics(Fingerprint_IDs, Plates_final, num_Fingerprint_per_plate, Reference_1, Reference_2):
    """
    Print statistics for every plate.

    Parameters
    ----------
    Fingerprint_IDs : array
        Unique Fingerprints.
    Plates_final : list of arrays
        Samples, References and Blanks distributed over num_plates plates.
    num_Fingerprint_per_plate : array
        Number of specific fingerprints on every plate.
    Reference_1 : str
        Label of Reference 1.
    Reference_2 : str
        Label of Reference 2.
    """

    print("-----------------------")  
    print("Summary:")  
    print("-----------------------")          
    for ind_p,plate in enumerate(Plates_final):
        print("Plate "+str(ind_p+1))
        print("#"+Reference_1+": %d" %len(np.where(plate==Reference_1)[0]))
        print("#"+Reference_2+": %d" %len(np.where(plate==Reference_2)[0]))
        print("#Blank: %d" %len(np.where(plate=="Blank")[0]))
        for index, fingerprint in enumerate(Fingerprint_IDs):
            print("#"+str(fingerprint)+":%d" %num_Fingerprint_per_plate[index,ind_p])   
        print("-----------------------")
    return

def get_Dictionaries():
    """
    Get dictionaries to translate between numbers and letters

    Returns
    -------
    dict_num2alph : dict
        Dictionary to translate numbers to letters.
    dict_alph2num : dict
        Dictionary to translate letters to numbers.
    """
    
    dict_num2alph={}
    dict_alph2num={}
    for num, alph in zip(range(26), string.ascii_uppercase):
        dict_num2alph.update({num:alph})
        dict_alph2num.update({alph:num})        
    return dict_num2alph, dict_alph2num

def read_plates(imp_columns, in_file="Out.txt", rows="Analysis.Row", columns="Analysis.Column", plate="Analysis.Plate"):
    """
    Read input file as written by generate_Output().

    Parameters
    ----------
    imp_columns : list of str
        Features used for Fingerprint determination.
    in_file : str, optional
        Input data (tab.separated) with headers. The default is "Out.txt".
    rows : str, optional
        Name of the column containing the row number. The default is "Analysis.Row".
    columns : str, optional
        Name of the column containing the column letter. The default is "Analysis.Column".
    plate : str, optional
        Name of the column containing the plate index. The default is "Analysis.Plate".

    Returns
    -------
    data : pandas.DataFrame
        Read-in data.
    num_rows : int
        Number of rows per plate.
    num_columns : int
        Number of columns per plate.  
    num_plates : int
        Number of plates.
    Fingerprint_IDs : array
        Unique Fingerprints.
    """
    
    dict_num2alph, dict_alph2num = get_Dictionaries()
    with open(in_file,"r") as file:
        data_raw = file.read().splitlines()
        data_raw = [line.replace("\"","") for line in data_raw]
        
    data_split = np.asarray([line.split(sep="\t") for line in data_raw])
    data = pd.DataFrame(data_split[1:] ,columns=data_split[0])
    
    num_rows = dict_alph2num[np.max(data[rows].values)]+1
    num_columns = np.max(np.int_(data[columns].values))
    num_plates = np.max(np.int_(data[plate].values))
    
    Fingerprints = []
    for index,row in data.iterrows(): 
        Fingerprints.append(row[imp_columns])
    
    # Generate Fingerprint
    Fingerprints = np.asarray(Fingerprints,dtype="U25")
    Fingerprint_IDs = np.unique(Fingerprints,axis=0)
    return data, num_rows, num_columns, num_plates, Fingerprint_IDs

def plot_PlateLayout(data, ref_column, imp_columns, Fingerprint_IDs, num_rows, num_columns, num_plates):
    """
    Plot the plate layout

    Parameters
    ----------
    data : pandas.DataFrame
        Data to be plotted.
    ref_column : str
        Label column.
    imp_columns : list of str
        Features used for Fingerprint determination.
    Fingerprint_IDs : array
        Unique Fingerprints.
    num_rows : int
        Number of rows per plate.
    num_columns : int
        Number of columns per plate.  
    num_plates : int
        Number of plates.
    """
    
    dict_num2alph, dict_alph2num = get_Dictionaries()
    Plates=np.zeros((num_rows,num_columns,num_plates),dtype="U25")
    Plates_color=np.zeros((num_rows,num_columns,num_plates))
    for index,row in data.iterrows():
        Plates[dict_alph2num[row["Analysis.Row"]],int(row["Analysis.Column"])-1,int(row["Analysis.Plate"])-1] = row[ref_column]
        for ID,fingerprint in enumerate(Fingerprint_IDs):
            if all(row[imp_columns]==fingerprint):
                    Plates_color[dict_alph2num[row["Analysis.Row"]],int(row["Analysis.Column"])-1,int(row["Analysis.Plate"])-1] = ID 
    
    colors_array=Color_Mix.colors(len(Fingerprint_IDs))
    cmap = colors.ListedColormap(colors_array)
    
    c_min = np.min(Plates_color)
    c_max = np.max(Plates_color)
    
    fs = 15
    
    for k in range(0,np.shape(Plates)[2],2):
        fig = plt.figure(constrained_layout=True)
        gs = fig.add_gridspec(2,1)
        fig.set_size_inches(20,16)
        rows = np.arange(num_rows-1,-1,-1)
        columns = np.arange(0,num_columns)
        
        ax1 = fig.add_subplot(gs[0,0])
        
        if k+1<np.shape(Plates)[2]:
            ax2 = fig.add_subplot(gs[1,0])
        
        for i, row in enumerate(rows):
            for j, column in enumerate(columns):
                c = Plates[i,j,k]
                ax1.text(column, row,c,va="center",ha="center")
        
        ax1.imshow(Plates_color[::-1,:,k],cmap=cmap,aspect=.5, vmin=c_min, vmax=c_max)
        
        ax1.set_xticks(columns)
        ax1.set_xticks(columns-.5, minor=True)
        ax1.set_xticklabels(np.int_(columns)+1,fontsize=fs)
        ax1.set_yticks(rows)
        ax1.set_yticks(rows-.5, minor=True)
        ax1.set_yticklabels([dict_num2alph[item] for item in rows][::-1],fontsize=fs)
        ax1.set_xlim(-.5,num_columns-.5)
        ax1.set_ylim(-.5,num_rows-.5)
        ax1.grid(which="minor",linestyle="--")  
        ax1.set_title("Plate "+str(k+1),fontsize=fs)
        
        if k+1<np.shape(Plates)[2]:
            for i, row in enumerate(rows):
                for j, column in enumerate(columns):
                    c = Plates[i,j,k+1]
                    ax2.text(column, row,c,va="center",ha="center")
            
            ax2.imshow(Plates_color[::-1,:,k+1],cmap=cmap,aspect=.5, vmin=c_min, vmax=c_max)
            
            ax2.set_xticks(columns)
            ax2.set_xticks(columns-.5, minor=True)
            ax2.set_xticklabels(np.int_(columns)+1,fontsize=fs)
            ax2.set_yticks(rows)
            ax2.set_yticks(rows-.5, minor=True)
            ax2.set_yticklabels([dict_num2alph[item] for item in rows][::-1],fontsize=fs)
            ax2.set_xlim(-.5,num_columns-.5)
            ax2.set_ylim(-.5,num_rows-.5)
            ax2.grid(which="minor",linestyle="--") 
            ax2.set_title("Plate "+str(k+2),fontsize=fs)
        plt.savefig("Plate"+str(k+1)+".png")  
    return