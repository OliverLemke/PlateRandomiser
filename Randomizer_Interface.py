# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Randomizer.ui'
#
# Created by: PyQt5 UI code generator 5.9.2

from PyQt5 import QtCore, QtGui, QtWidgets
import numpy as np
import pandas as pd
from sklearn.utils import shuffle
import string
import Color_Mix
import matplotlib.pyplot as plt
from matplotlib import colors
import os
from sklearn.metrics.pairwise import manhattan_distances

### To Do
# Add forbidden/fixed wells

# Name
# Icon

# Select plate Layout
def get_Plate_layout(Plate_layout, num_Blanks, num_References=None, num_Ref_1=None, num_Ref_2=None):
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
def get_data_dicts(data, ref_columns, imp_columns):
    """
    Generates Fingerprint dictionary

    Parameters
    ----------
    data : pandas.DataFrame
        Input data.
    ref_columns : list of str
        Identifier columns for the sample.
    imp_columns : list of str
        Features used for Fingerprint determination.

    Returns
    -------
    data_dict : dict
        Dictionary containing the data with ref_columns[0] as key.
    Fingerprints : array
        Fingerprints for all samples (row-wise).

    """
    data_dict = {}
    Fingerprints = []
    for index,row in data.iterrows():
        data_dict.update({row[ref_columns[0]]:row})
        Fingerprints.append(row[imp_columns])
    return data_dict, Fingerprints

def get_fixed_stats(data, ref_column, label_column, cut):
    """
    Get Dictionaries for relationships

    Parameters
    ----------
    data : pandas.DataFrame
        Input data.
    ref_columns : list of str
        Identifier columns for the sample.
    label_column : str
        Column used as fixed references.
    cut : int
        Cutoff for pred-istirbuted samples.

    Returns
    -------
    dict_fix_to_label : dict
        Dictionary for fixed references to label.
    dict_label_to_fix : dict
        Dictionary for label to references column.
    dict_fix_to_fingerprint : dict
        Dictionary for references column to fingerprint.
    labels_cut : list of str
        Samples to be pre-distributed.
    """
    
    # B = dict_fix_to_label
    dict_fix_to_label = {}
    for ID in data[label_column].values:
        dict_fix_to_label.update({ID:data[data[label_column]==ID][ref_column].values})
    # C = dict_label_to_fix
    dict_label_to_fix = {}
    for ID in data[label_column].values:
        for item in data[data[label_column]==ID][ref_column].values:
            dict_label_to_fix.update({item:ID})
      
    # A = dict_fix_to_fingerprint
    dict_fix_to_fingerprint = {}
    # D = labels_cut
    labels_cut = []
    
    for key in np.unique(data[label_column].values):
        dict_fix_to_fingerprint.update({key:len(data[data[label_column]==key][label_column].values)})
        if len(data[data[label_column]==key][label_column].values) >= cut:
            labels_cut.append(key)
    return dict_fix_to_label, dict_label_to_fix, dict_fix_to_fingerprint, labels_cut

# Read data: add tsv, csv and xlsx
def get_data_dict(data, ref_columns, imp_columns):
    """
    Generates Fingerprint dictionary
    Parameters
    ----------
    data : pandas.DataFrame
        Input data.
    ref_columns : list of str
        Identifier columns for the sample.
    imp_columns : list of str
        Features used for Fingerprint determination.
    Returns
    -------
    data_dict : dict
        Dictionary containing the data with ref_columns[0] as key.
    Fingerprints : array
        Fingerprints for all samples (row-wise).
    """
    data_dict = {}
    Fingerprints = []
    for index,row in data.iterrows():
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
    num_Fingerprint: array
        Number of Occurence for every Fingerprint.
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
    num_Fingerprint = np.asarray([len(item) for item in Fingerprints_list])
    return Fingerprint_IDs, Fingerprints_list, num_Fingerprint

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
        print("Number of Plates set to %d" %num_plates)
    
    ## Determine how many per plate (floor)
    num_Sample_per_plate = np.int_(np.floor(num_Fingerprint/num_plates))
    
    ## Determine how many are too large (mod)
    num_Samples_overlap = np.int_(np.mod(num_Fingerprint,num_plates))
    return num_Fingerprint, num_plates, num_Sample_per_plate, num_Samples_overlap

# Distribute Samples across the plates
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
        plate = np.asarray([], dtype='<U25')
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

# Distribute References across the plates
def distribute_References(Plates_overlap, Fingerprints_list, Reference_1, Reference_2, num_wells, num_Blanks, num_Ref_1=6, num_Ref_2=6, percentage_Ref_1=None, dict_fix_to_label=None):
    """
    Add references and blanks to the plate

    Parameters
    ----------
    Plates_overlap : list of arrays
        Samples distributed num_plates over plates.
    Fingerprints_list : list of str
        Samples assigned to each Fingerprint.
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
        plate = shuffle(np.concatenate((plate,Ref_1,Ref_2,Blanks)))
        score = len(plate)**2
        for i in range(250):
            plate = shuffle(plate)
            score_eval = evaluate_plate_layout(plate, Fingerprints_list, Reference_1, Reference_2)
            if dict_fix_to_label:
                score_label = evaluate_plate_layout_fixed(plate,dict_fix_to_label)
                score_eval = (score_eval+score_label)/2
            if score_eval <= score:
                score = score_eval
                print(i)
                print(score)
                plate_final = plate
        Plates_final.append(plate_final)
    return Plates_final

def evaluate_plate_layout_fixed(plate,dict_fix_to_label):
    score = 0
    if len(plate)==96:
        num_columns = 12
    elif len(plate)==384:
        num_columns = 24
    for key in dict_fix_to_label:
        coord = np.asarray([[int(ind/num_columns),np.mod(ind,num_columns)] for ind,item in enumerate(plate) if item.astype("U25") in dict_fix_to_label[key].astype("U25")])
        try:
            score += (np.sum(np.exp(-1/len(coord)*manhattan_distances(coord)))-len(coord))/2
        except:
            pass
    return score

def get_Fingerprints_list_full_fixed(Fingerprints_list, dict_fix_to_label):
    Fingerprints_list_full = [np.asarray([el for key in item for el in dict_fix_to_label[key]]) for item in Fingerprints_list]
    return Fingerprints_list_full

def evaluate_plate_layout(plate, Fingerprints_list, Reference_1, Reference_2):
    score = 0
    if len(plate)==96:
        num_columns = 12
    elif len(plate)==384:
        num_columns = 24
    for el in Fingerprints_list:
        coord = np.asarray([[int(ind/num_columns),np.mod(ind,num_columns)] for ind,item in enumerate(plate) if item.astype("U25") in el.astype("U25")])
        try:
            score += (np.sum(np.exp(-1*manhattan_distances(coord)))-len(coord))/2
        except:
            pass
    coord = np.asarray([[int(ind/num_columns),np.mod(ind,num_columns)] for ind,item in enumerate(plate) if item == Reference_1])
    try:
        score += (np.sum(np.exp(-1*manhattan_distances(coord)))-len(coord))/2
    except:
        pass
    coord = np.asarray([[int(ind/num_columns),np.mod(ind,num_columns)] for ind,item in enumerate(plate) if item == Reference_2])
    try:
        score += (np.sum(np.exp(-1*manhattan_distances(coord)))-len(coord))/2
    except:
        pass
    return score
    ##### negative exp ### exp(-D)

# Wirte Output-file
def generate_Output(data_dict, Plates_final, ref_columns, imp_columns, num_columns, remaining_columns, Reference_1, Reference_2, out_file="Out.txt", label_column = None):
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
    output_columns : list of str
        Additional columns to be added to output.
    out_file : str, optional
        Name of the output file. The default is "Out.txt".
    label_column : list of str, optional
        List containing the column used for force samples onto one plate. The default is None.
    """

    dict_alph, dict_num = get_Dictionaries()
    
    c_Ref_1 = 1
    c_Ref_2 = 1
    c_Blank = 1
    
    if out_file[-4:]==".csv":
        sep=","
    else:
        sep="\t"
    
    with open(out_file,"w") as file:
        file.write("Analysis.Plate"+sep+"Analysis.Row"+sep+"Analysis.Column"+sep+"Sample.Name"+sep+"Peptide.Concentration"+sep)
        for ref_column in ref_columns:
            file.write("%s" %ref_column)
            file.write(sep)
        for column in imp_columns:
            file.write("%s" %column)
            file.write(sep)
        if label_column:
            file.write("%s" %label_column[0])
            file.write(sep)
            remaining_columns = [column for column in remaining_columns if column not in label_column]
        
        remaining_columns = [column for column in remaining_columns if (column not in ref_columns) & (column not in imp_columns)]
        
        if remaining_columns==[]:        
            file.write("\n")
        else:
            for ind,column in enumerate(remaining_columns):
                if ind != len(remaining_columns)-1:
                    file.write("%s" %column)
                    file.write(sep)
                else:
                    file.write("%s\n" %column)    
        for index,plate in enumerate(Plates_final):
            for index2,el in enumerate(plate):
                file.write(str(index+1)+sep)
                file.write("%s" %dict_alph[np.floor(index2/num_columns)])
                file.write(sep)
                file.write(str(np.mod(index2,num_columns)+1)+sep)
                if (el != Reference_1) & (el != Reference_2) & (el != "Blank"):
                    try:
                        check = data_dict[el]
                    except:
                        check = data_dict[int(float(el))]
                        el = int(float(el))
                if el == Reference_1:
                    file.write(Reference_1+"_{:03}".format(c_Ref_1))
                elif el == Reference_2:
                    file.write(Reference_2+"_{:03}".format(c_Ref_2))
                elif el == "Blank":
                    file.write("Blank_"+"{:03}".format(c_Blank))
                else:
                    file.write(str(el))
                file.write(sep)
                file.write("1"+sep)
                if el == Reference_1:
                    for i in range(len(ref_columns)):
                        file.write(Reference_1+"_{:03}".format(c_Ref_1)+sep)
                    for i in range(len(imp_columns)):
                        file.write("QC_1")
                        file.write(sep)
                    if label_column:
                        file.write("QC_1")
                        file.write(sep)
                    c_Ref_1+=1
                elif el == Reference_2:
                    for i in range(len(ref_columns)):
                        file.write(Reference_2+"_{:03}".format(c_Ref_2)+sep)
                    for i in range(len(imp_columns)):
                        file.write("QC_2")
                        file.write(sep)
                    if label_column:
                        file.write("QC_2")
                        file.write(sep)
                    c_Ref_2+=1
                elif el == "Blank":
                    for i in range(len(ref_columns)):
                        file.write("Blank_"+"{:03}".format(c_Blank)+sep)
                    for i in range(len(imp_columns)):
                        file.write("Blank")
                        file.write(sep)
                    if label_column:
                        file.write("Blank")
                        file.write(sep)
                    c_Blank+=1
                else:
                    for ref_column in ref_columns:
                        file.write("%s" %data_dict[el][ref_column])
                        file.write(sep)
                    for column in imp_columns:
                        file.write("%s" %data_dict[el][column])
                        file.write(sep)
                    if label_column:
                        file.write("%s" %data_dict[el][label_column[0]])
                        file.write(sep)
                if remaining_columns==[]:        
                    file.write("\n")
                else:
                    for ind,column in enumerate(remaining_columns):
                        if ind != len(remaining_columns)-1:
                            if (el != Reference_1) & (el != Reference_2) & (el != "Blank"):
                                file.write("%s" %data_dict[el][column])
                                file.write(sep)
                            else:
                                file.write("NA")
                                file.write(sep)
                        else:
                            if (el != Reference_1) & (el != Reference_2) & (el != "Blank"):
                                file.write("%s\n" %data_dict[el][column])
                            else:
                                file.write("NA\n")
    return

# Get Plate statistics
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
            try:
                check = data_dict[el]
            except:
                try:
                    check = data_dict[int(float(el))]
                    el = int(float(el))
                except:
                    pass
            for index,fingerprint in enumerate(Fingerprint_IDs):
                try:
                    if all(data_dict[el][imp_columns].values.astype("U25")==fingerprint):
                        num_Fingerprint_per_plate[index,ind_p]+=1
                except:
                    pass
    return num_Fingerprint_per_plate

# Print Plate statistics
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

# Get Dictionaries Alph <-> Num
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

# Read Outputfile for plotting
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
    
    if in_file[:-4]==".csv":
        sep=","
    else:
        sep="\t"
    
    dict_num2alph, dict_alph2num = get_Dictionaries()
    with open(in_file,"r") as file:
        data_raw = file.read().splitlines()
        data_raw = [line.replace("\"","") for line in data_raw]
        
    data_split = np.asarray([line.split(sep=sep) for line in data_raw])
    data = pd.DataFrame(data_split[1:] ,columns=data_split[0])
    
    Fingerprints = []
    for index,row in data.iterrows(): 
        Fingerprints.append(row[imp_columns])
    
    # Generate Fingerprint
    Fingerprints = np.asarray(Fingerprints,dtype="U25")
    Fingerprint_IDs = np.unique(Fingerprints,axis=0)
    return data, Fingerprint_IDs

# Plot Plate-layout
def plot_PlateLayout(data, ref_column, imp_columns, Fingerprint_IDs, num_rows, num_columns, num_plates, out_file, label_column=None):
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
    out_file : str
        Name of the output file.
    label_column : str, optional
        Name of the column that was used to force samples onto the same plate. The default is None.
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
    
    for k in range(0,num_plates,2):
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
        plt.savefig(out_file+str(k+1)+"_and_"+str(k+2)+".png")  
        
    for k in range(0,num_plates):
        fig,ax = plt.subplots()
        fig.set_size_inches(20,8)
        rows = np.arange(num_rows-1,-1,-1)
        columns = np.arange(0,num_columns)
        
        for i, row in enumerate(rows):
            for j, column in enumerate(columns):
                c = Plates[i,j,k]
                ax.text(column, row,c,va="center",ha="center")
        
        ax.imshow(Plates_color[::-1,:,k],cmap=cmap,aspect=.5, vmin=c_min, vmax=c_max)
        
        ax.set_xticks(columns)
        ax.set_xticks(columns-.5, minor=True)
        ax.set_xticklabels(np.int_(columns)+1,fontsize=fs)
        ax.set_yticks(rows)
        ax.set_yticks(rows-.5, minor=True)
        ax.set_yticklabels([dict_num2alph[item] for item in rows][::-1],fontsize=fs)
        ax.set_xlim(-.5,num_columns-.5)
        ax.set_ylim(-.5,num_rows-.5)
        ax.grid(which="minor",linestyle="--")  
        ax.set_title("Plate "+str(k+1),fontsize=fs)
        plt.savefig(out_file+str(k+1)+".png")

    if label_column:
        Plates=np.zeros((num_rows,num_columns,num_plates),dtype="U25")
        for index,row in data.iterrows():
            Plates[dict_alph2num[row["Analysis.Row"]],int(row["Analysis.Column"])-1,int(row["Analysis.Plate"])-1] = row[label_column]

        for k in range(0,num_plates):
            fig,ax = plt.subplots()
            fig.set_size_inches(20,8)
            rows = np.arange(num_rows-1,-1,-1)
            columns = np.arange(0,num_columns)
            
            for i, row in enumerate(rows):
                for j, column in enumerate(columns):
                    c = Plates[i,j,k]
                    ax.text(column, row,c,va="center",ha="center")
            
            ax.imshow(Plates_color[::-1,:,k],cmap=cmap,aspect=.5, vmin=c_min, vmax=c_max)
            
            ax.set_xticks(columns)
            ax.set_xticks(columns-.5, minor=True)
            ax.set_xticklabels(np.int_(columns)+1,fontsize=fs)
            ax.set_yticks(rows)
            ax.set_yticks(rows-.5, minor=True)
            ax.set_yticklabels([dict_num2alph[item] for item in rows][::-1],fontsize=fs)
            ax.set_xlim(-.5,num_columns-.5)
            ax.set_ylim(-.5,num_rows-.5)
            ax.grid(which="minor",linestyle="--")  
            ax.set_title("Plate "+str(k+1),fontsize=fs)
            plt.savefig(out_file+str(k+1)+"_label.png")  
        return

# Get Number of Fingerprints for forced samples
def get_num_Fingerprint(Fingerprints_list, num_wells_to_fill, dict_fix_to_fingerprint):
    """
    Get Number of Fingerprints with forcing samples onto one plate.

    Parameters
    ----------
    Fingerprints_list : list of str
        Samples assigned to each Fingerprint.
    num_wells_to_fill : int
        Maximal number of wells to be filled per plate.
    dict_fix_to_fingerprint : dict
        Dictionary for references column to fingerprint.

    Returns
    -------
    num_Fingerprint : array
        Number of samples per fingerprint
    num_plates : int
        Final number of plates.
    """
    
    num_Fingerprint = np.asarray([np.sum(np.asarray([dict_fix_to_fingerprint[key] for key in item])) for item in Fingerprints_list])
    num_plates = int(np.ceil(np.sum(num_Fingerprint)/num_wells_to_fill))
    return num_Fingerprint, num_plates

# Distribute samples with forcing samples onto the same plate
def distribute_samples_fixed(Fingerprints_list, data_dict, num_plates, dict_fix_to_label, labels_cut, num_wells, num_References, num_Blanks):
    """
    Distribute samples forcing samples of same identifier onto the same plate.

    Parameters
    ----------
    Fingerprints_list : list of str
        Samples assigned to each Fingerprint.
    data_dict : dict
        Dictionary containing the data with ref_columns[0] as key.
    num_plates : int
        Final number of plates.
    dict_fix_to_label : dict
        Dictionary for fixed references to label.
    labels_cut : list of str
        Samples to be pre-distributed.
    num_wells : int
        Number of wells per plate.
    num_References : int
        Number of References.
    num_Blanks : int
        Number of blanks to be included per plate.

    Returns
    -------
    Plates_full : list of arrays
        Samples distributed over num_plates plates.
    Plates : list of arrays
        Samples distributed over num_plates plates. Identifier used for fixing data onto plates.
    """

    Repeat = True
    while Repeat:        
        Fingerprints_list_new = [shuffle([item for item in fingerprint if item not in labels_cut]) for fingerprint in Fingerprints_list]           
        ## Determine how many per plate (floor)
        num_Sample_per_plate = np.int_(np.floor(np.asarray([len(item) for item in Fingerprints_list_new])/num_plates))
        ## Determine how many are too large (mod)
        num_Samples_overlap = np.int_(np.mod(np.asarray([len(item) for item in Fingerprints_list_new]),num_plates))
        
        Plates = []
        labels_cut = shuffle(np.asarray(labels_cut))
        
        for i in range(num_plates):
            plate = np.asarray(labels_cut[i::num_plates])            
            for ind,fingerprint in enumerate(Fingerprints_list_new):
                ## Cut at pre-determined cut
                plate = np.concatenate((plate,fingerprint[i*num_Sample_per_plate[ind]:(i+1)*num_Sample_per_plate[ind]]))
            Plates.append(plate)
        
        Distributed_samples = [item for plate in Plates for item in plate]
        All_samples = [key for key in dict_fix_to_label]
            
        ## Distribute remaining data points        
        Overlap = [item for item in All_samples if item not in Distributed_samples]        
        Overlap = shuffle(Overlap)
        Plates_step = [[el for item in plate for el in dict_fix_to_label[item]] for plate in Plates]
        
        if all(np.asarray([len(plate) for plate in Plates_step]) <= num_wells-num_References-num_Blanks):
            c=0
            b=0
            while c<len(Overlap):
                Plates_step = [[el for item in plate for el in dict_fix_to_label[item]] for plate in Plates]
                for ind,plate in enumerate(Plates):
                    if len(Plates_step[ind]) + len(dict_fix_to_label[Overlap[c]]) <= num_wells-num_Blanks-num_References:
                        Plates[ind]=(np.concatenate((plate,[Overlap[c]])))
                        c+=1
                    b+=1
                    if c == len(Overlap):
                        break
                if b>=200:
                    print("failed")
                    break
            if c == len(Overlap):    
                Plates_full = [[el for item in plate for el in dict_fix_to_label[item]] for plate in Plates]
            else:
                Plates_full = [np.zeros(100)]
        else: 
            Plates_full = [np.zeros(100)]
        
        if any(np.asarray([len(plate) for plate in Plates_full]) > num_wells-num_References-num_Blanks):
            #print(np.asarray([len(plate) for plate in Plates_full]))
            print("try again")
        else:
            Repeat = False
            
    Check = np.asarray([item for plate in Plates for item in plate])
    print("Unique IDs provided: %d" %len(data_dict))
    print("Unique IDs detected on plates: %d" %len(np.unique(Check)))
    print("Total IDs on plates: %d" %(np.sum(np.asarray([len(item) for item in Plates]))))
    print("Unique IDs on plates: %d" %(len(np.unique(np.asarray([el for item in Plates for el in item])))))
    return Plates_full, Plates

# Dictionary for data if samples were forced
def get_data_dict_fixed(data, ref_columns):
    """
    Generate dictionary for forced samples.

    Parameters
    ----------
    data : pandas.DataFrame
        Input data.
    ref_columns : list of str
        Identifier columns for the sample.

    Returns
    -------
    data_dict : dict
        Dictionary containing the data with ref_columns[0] as key.
    """
    
    data_dict_fixed = {}
    for index,row in data.iterrows():
        data_dict_fixed.update({row[ref_columns[0]]:row})
    return data_dict_fixed

#%%

class Stream(QtCore.QObject):
    newText = QtCore.pyqtSignal(str)

    def write(self, text):
        self.newText.emit(str(text))
    
    def flush(self):
        pass

class Ui_Form(QtWidgets.QWidget):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(1273, 689)
        self.groupBox_Input = QtWidgets.QGroupBox(Form)
        self.groupBox_Input.setGeometry(QtCore.QRect(10, 10, 311, 671))
        self.groupBox_Input.setObjectName("groupBox_Input")
        self.horizontalLayoutWidget = QtWidgets.QWidget(self.groupBox_Input)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(10, 20, 291, 41))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout_Input = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout_Input.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_Input.setObjectName("horizontalLayout_Input")
        self.lineEdit_Input = QtWidgets.QLineEdit(self.horizontalLayoutWidget)
        self.lineEdit_Input.setObjectName("lineEdit_Input")
        self.horizontalLayout_Input.addWidget(self.lineEdit_Input)
        self.pushButton_Input = QtWidgets.QPushButton(self.horizontalLayoutWidget)
        self.pushButton_Input.setObjectName("pushButton_Input")
        self.horizontalLayout_Input.addWidget(self.pushButton_Input)
        self.pushButton_Recompute = QtWidgets.QPushButton(self.groupBox_Input)
        self.pushButton_Recompute.setGeometry(QtCore.QRect(190, 60, 113, 32))
        self.pushButton_Recompute.setObjectName("pushButton_Recompute")
        self.scrollArea_Sheets = QtWidgets.QScrollArea(self.groupBox_Input)
        self.scrollArea_Sheets.setGeometry(QtCore.QRect(10, 90, 291, 121))
        self.scrollArea_Sheets.setWidgetResizable(True)
        self.scrollArea_Sheets.setObjectName("scrollArea_Sheets")
        self.scrollAreaWidgetContents_Sheets = QtWidgets.QWidget() #QButtonGroup(self, exclusive=False)
        self.scrollAreaWidgetContents_Sheets.setGeometry(QtCore.QRect(0, 0, 289, 119))
        self.scrollAreaWidgetContents_Sheets.setObjectName("scrollAreaWidgetContents_Sheets")
        self.scrollArea_Sheets.setWidget(self.scrollAreaWidgetContents_Sheets)
        self.label_inp_columns = QtWidgets.QLabel(self.groupBox_Input)
        self.label_inp_columns.setGeometry(QtCore.QRect(10, 220, 291, 16))
        self.label_inp_columns.setObjectName("label_inp_columns")
        self.scrollArea_imp_columns = QtWidgets.QScrollArea(self.groupBox_Input)
        self.scrollArea_imp_columns.setGeometry(QtCore.QRect(10, 240, 291, 121))
        self.scrollArea_imp_columns.setWidgetResizable(True)
        self.scrollArea_imp_columns.setObjectName("scrollArea_imp_columns")
        self.scrollAreaWidgetContents_imp_columns = QtWidgets.QWidget()
        self.scrollAreaWidgetContents_imp_columns.setGeometry(QtCore.QRect(0, 0, 289, 119))
        self.scrollAreaWidgetContents_imp_columns.setObjectName("scrollAreaWidgetContents_imp_columns")
        self.scrollArea_imp_columns.setWidget(self.scrollAreaWidgetContents_imp_columns)
        self.label_output_columns = QtWidgets.QLabel(self.groupBox_Input)
        self.label_output_columns.setGeometry(QtCore.QRect(10, 520, 291, 16))
        self.label_output_columns.setObjectName("label_output_columns")
        self.label_ref_columns = QtWidgets.QLabel(self.groupBox_Input)
        self.label_ref_columns.setGeometry(QtCore.QRect(10, 370, 291, 16))
        self.label_ref_columns.setObjectName("label_ref_columns")
        self.scrollArea_output_columns = QtWidgets.QScrollArea(self.groupBox_Input)
        self.scrollArea_output_columns.setGeometry(QtCore.QRect(10, 540, 291, 121))
        self.scrollArea_output_columns.setWidgetResizable(True)
        self.scrollArea_output_columns.setObjectName("scrollArea_output_columns")
        self.scrollAreaWidgetContents_output_columns = QtWidgets.QWidget()
        self.scrollAreaWidgetContents_output_columns.setGeometry(QtCore.QRect(0, 0, 289, 119))
        self.scrollAreaWidgetContents_output_columns.setObjectName("scrollAreaWidgetContents_output_columns")
        self.scrollArea_output_columns.setWidget(self.scrollAreaWidgetContents_output_columns)
        self.label_Sheets = QtWidgets.QLabel(self.groupBox_Input)
        self.label_Sheets.setGeometry(QtCore.QRect(10, 70, 171, 16))
        self.label_Sheets.setObjectName("label_Sheets")
        self.scrollArea_ref_columns = QtWidgets.QScrollArea(self.groupBox_Input)
        self.scrollArea_ref_columns.setGeometry(QtCore.QRect(10, 390, 291, 121))
        self.scrollArea_ref_columns.setWidgetResizable(True)
        self.scrollArea_ref_columns.setObjectName("scrollArea_ref_columns")
        self.scrollAreaWidgetContents_ref_columns = QtWidgets.QWidget()
        self.scrollAreaWidgetContents_ref_columns.setGeometry(QtCore.QRect(0, 0, 289, 119))
        self.scrollAreaWidgetContents_ref_columns.setObjectName("scrollAreaWidgetContents_ref_columns")
        self.scrollArea_ref_columns.setWidget(self.scrollAreaWidgetContents_ref_columns)
        
        self.groupBox_Fingerprints_Preview = QtWidgets.QGroupBox(Form)
        self.groupBox_Fingerprints_Preview.setGeometry(QtCore.QRect(330, 10, 301, 361))
        self.groupBox_Fingerprints_Preview.setObjectName("groupBox_Fingerprints_Preview")
        self.scrollArea_Fingerprints_Preview = QtWidgets.QScrollArea(self.groupBox_Fingerprints_Preview)
        self.scrollArea_Fingerprints_Preview.setGeometry(QtCore.QRect(10, 20, 281, 301))
        self.scrollArea_Fingerprints_Preview.setWidgetResizable(True)
        self.scrollArea_Fingerprints_Preview.setObjectName("scrollArea_Fingerprints_Preview")
        self.scrollAreaWidgetContents_Fingerprint_Preview = QtWidgets.QWidget()
        self.scrollAreaWidgetContents_Fingerprint_Preview.setGeometry(QtCore.QRect(0, 0, 279, 299))
        self.scrollAreaWidgetContents_Fingerprint_Preview.setObjectName("scrollAreaWidgetContents_Fingerprint_Preview")
        self.plainTextEdit_Fingerprints_Preview = QtWidgets.QPlainTextEdit(self.scrollAreaWidgetContents_Fingerprint_Preview)
        #self.plainTextEdit_Fingerprints_Preview.setEnabled(True)
        self.plainTextEdit_Fingerprints_Preview.setGeometry(QtCore.QRect(0, 0, 281, 301))
        self.plainTextEdit_Fingerprints_Preview.setReadOnly(True)
        self.plainTextEdit_Fingerprints_Preview.setPlainText("")
        self.plainTextEdit_Fingerprints_Preview.setObjectName("plainTextEdit_Fingerprints_Preview")
        self.scrollArea_Fingerprints_Preview.setWidget(self.scrollAreaWidgetContents_Fingerprint_Preview)
        self.pushButton_Fingerprints_Preview = QtWidgets.QPushButton(self.groupBox_Fingerprints_Preview)
        self.pushButton_Fingerprints_Preview.setGeometry(QtCore.QRect(12, 330, 281, 32))
        self.pushButton_Fingerprints_Preview.setObjectName("pushButton_Fingerprints_Preview")
        
        self.groupBox_Parameters = QtWidgets.QGroupBox(Form)
        self.groupBox_Parameters.setGeometry(QtCore.QRect(330, 370, 301, 311))
        self.groupBox_Parameters.setObjectName("groupBox_Parameters")
        self.gridLayoutWidget = QtWidgets.QWidget(self.groupBox_Parameters)
        self.gridLayoutWidget.setGeometry(QtCore.QRect(10, 20, 281, 121))
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridlayout = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridlayout.setContentsMargins(0, 0, 0, 0)
        self.gridlayout.setObjectName("gridlayout")
        self.label_Layout = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_Layout.setObjectName("label_Layout")
        self.gridlayout.addWidget(self.label_Layout, 0, 0, 1, 1)
        self.comboBox_Layout = QtWidgets.QComboBox(self.gridLayoutWidget)
        self.comboBox_Layout.setObjectName("comboBox_Layout")
        self.comboBox_Layout.addItems(["96","384"])
        self.gridlayout.addWidget(self.comboBox_Layout, 0, 1, 1, 1)
        self.label_Blanks = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_Blanks.setObjectName("label_Blanks")
        self.gridlayout.addWidget(self.label_Blanks, 1, 0, 1, 1)
        self.spinBox_Blanks = QtWidgets.QSpinBox(self.gridLayoutWidget)
        self.spinBox_Blanks.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.spinBox_Blanks.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.spinBox_Blanks.setProperty("value", 1)
        self.spinBox_Blanks.setObjectName("spinBox_Blanks")
        self.gridlayout.addWidget(self.spinBox_Blanks, 1, 1, 1, 1)
        self.label_Ref1 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_Ref1.setObjectName("label_Ref1")
        self.gridlayout.addWidget(self.label_Ref1, 2, 0, 1, 1)
        self.lineEdit_Ref1 = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_Ref1.setObjectName("lineEdit_Ref1")
        self.lineEdit_Ref1.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.gridlayout.addWidget(self.lineEdit_Ref1, 2, 1, 1, 1)
        self.label_Ref2 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_Ref2.setObjectName("label_Ref2")
        self.gridlayout.addWidget(self.label_Ref2, 3, 0, 1, 1)
        self.lineEdit_Ref2 = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.lineEdit_Ref2.setObjectName("lineEdit_Ref2")
        self.lineEdit_Ref2.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.gridlayout.addWidget(self.lineEdit_Ref2, 3, 1, 1, 1)
        
        self.tabWidget = QtWidgets.QTabWidget(self.groupBox_Parameters)
        self.tabWidget.setGeometry(QtCore.QRect(10, 150, 281, 101))
        self.tabWidget.setObjectName("tabWidget")
        self.tab_Fixed = QtWidgets.QWidget()
        self.tab_Fixed.setObjectName("tab_Fixed")
        self.gridLayoutWidget_2 = QtWidgets.QWidget(self.tab_Fixed)
        self.gridLayoutWidget_2.setGeometry(QtCore.QRect(0, 0, 271, 71))
        self.gridLayoutWidget_2.setObjectName("gridLayoutWidget_2")
        self.gridLayout_Fixed = QtWidgets.QGridLayout(self.gridLayoutWidget_2)
        self.gridLayout_Fixed.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_Fixed.setObjectName("gridLayout_Fixed")
        self.label_num_Ref1 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_num_Ref1.setObjectName("label_num_Ref1")
        self.gridLayout_Fixed.addWidget(self.label_num_Ref1, 0, 0, 1, 1)
        self.spinBox_num_Ref1 = QtWidgets.QSpinBox(self.gridLayoutWidget_2)
        self.spinBox_num_Ref1.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.spinBox_num_Ref1.setMaximum(97)
        self.spinBox_num_Ref1.setObjectName("spinBox_num_Ref1")
        self.gridLayout_Fixed.addWidget(self.spinBox_num_Ref1, 0, 1, 1, 1)
        self.label_num_Ref2 = QtWidgets.QLabel(self.gridLayoutWidget_2)
        self.label_num_Ref2.setObjectName("label_num_Ref2")
        self.gridLayout_Fixed.addWidget(self.label_num_Ref2, 1, 0, 1, 1)
        self.spinBox_num_Ref2 = QtWidgets.QSpinBox(self.gridLayoutWidget_2)
        self.spinBox_num_Ref2.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.spinBox_num_Ref2.setObjectName("spinBox_num_Ref2")
        self.gridLayout_Fixed.addWidget(self.spinBox_num_Ref2, 1, 1, 1, 1)
        self.tabWidget.addTab(self.tab_Fixed, "")
        self.tab_Relative = QtWidgets.QWidget()
        self.tab_Relative.setObjectName("tab_Relative")
        self.gridLayoutWidget_3 = QtWidgets.QWidget(self.tab_Relative)
        self.gridLayoutWidget_3.setGeometry(QtCore.QRect(0, 0, 271, 71))
        self.gridLayoutWidget_3.setObjectName("gridLayoutWidget_3")
        self.gridLayout_Relative = QtWidgets.QGridLayout(self.gridLayoutWidget_3)
        self.gridLayout_Relative.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_Relative.setObjectName("gridLayout_Relative")
        self.spinBox_min_Ref = QtWidgets.QSpinBox(self.gridLayoutWidget_3)
        self.spinBox_min_Ref.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.spinBox_min_Ref.setObjectName("spinBox_min_Ref")
        self.gridLayout_Relative.addWidget(self.spinBox_min_Ref, 0, 1, 1, 1)
        self.label_min_Ref = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_min_Ref.setObjectName("label_min_Ref")
        self.gridLayout_Relative.addWidget(self.label_min_Ref, 0, 0, 1, 1)
        self.label_perc_Ref1 = QtWidgets.QLabel(self.gridLayoutWidget_3)
        self.label_perc_Ref1.setObjectName("label_perc_Ref1")
        self.gridLayout_Relative.addWidget(self.label_perc_Ref1, 1, 0, 1, 1)
        self.doubleSpinBox_perc_Ref1 = QtWidgets.QDoubleSpinBox(self.gridLayoutWidget_3)
        self.doubleSpinBox_perc_Ref1.setMaximum(1.0)
        self.doubleSpinBox_perc_Ref1.setSingleStep(0.01)
        self.doubleSpinBox_perc_Ref1.setProperty("value", 1.0)
        self.doubleSpinBox_perc_Ref1.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.doubleSpinBox_perc_Ref1.setObjectName("doubleSpinBox_perc_Ref1")
        self.gridLayout_Relative.addWidget(self.doubleSpinBox_perc_Ref1, 1, 1, 1, 1)
        self.tabWidget.addTab(self.tab_Relative, "")
        
        self.horizontalLayoutWidget_Seed = QtWidgets.QWidget(self.groupBox_Parameters)
        self.horizontalLayoutWidget_Seed.setGeometry(QtCore.QRect(0, 280, 301, 24))
        self.horizontalLayoutWidget_Seed.setObjectName("horizontalLayoutWidget_Seed")
        self.horizontalLayout_Seed = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_Seed)
        self.horizontalLayout_Seed.setContentsMargins(10, 0, 10, 0)
        self.horizontalLayout_Seed.setObjectName("horizontalLayout_Seed")
        #self.horizontalSlider_Seed = QtWidgets.QSlider(self.horizontalLayoutWidget_Seed)
        #self.horizontalSlider_Seed.setMaximum(99999)
        #self.horizontalSlider_Seed.setSliderPosition(42)
        #self.horizontalSlider_Seed.setOrientation(QtCore.Qt.Horizontal)
        #self.horizontalSlider_Seed.setObjectName("horizontalSlider_Seed")
        #self.horizontalLayout_Seed.addWidget(self.horizontalSlider_Seed)
        self.lineEdit_Seed = QtWidgets.QLineEdit(self.horizontalLayoutWidget_Seed)
        self.lineEdit_Seed.setObjectName("lineEdit_Seed")
        self.lineEdit_Seed.setText("42")
        self.lineEdit_Seed.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.horizontalLayout_Seed.addWidget(self.lineEdit_Seed)
        self.label_Seed = QtWidgets.QLabel(self.groupBox_Parameters)
        self.label_Seed.setGeometry(QtCore.QRect(10, 260, 60, 16))
        self.label_Seed.setObjectName("label_Seed")
        self.checkbox_fix_column = QtWidgets.QCheckBox(self.horizontalLayoutWidget_Seed)
        self.checkbox_fix_column.setObjectName("checkbox_fix_column")
        self.horizontalLayout_Seed.addWidget(self.checkbox_fix_column)
        self.checkbox_optimize = QtWidgets.QCheckBox(self.horizontalLayoutWidget_Seed)
        self.checkbox_optimize.setObjectName("checkbox_optimize")
        self.horizontalLayout_Seed.addWidget(self.checkbox_optimize)
        self.checkbox_exclude = QtWidgets.QCheckBox(self.horizontalLayoutWidget_Seed)
        self.checkbox_exclude.setObjectName("checkbox_exclude")
        self.horizontalLayout_Seed.addWidget(self.checkbox_exclude)
        self.label_Settings = QtWidgets.QLabel(self.groupBox_Parameters)
        self.label_Settings.setGeometry(QtCore.QRect(131, 260, 60, 16))
        self.label_Settings.setObjectName("label_Settings")
        
        self.groupBox_Output = QtWidgets.QGroupBox(Form)
        self.groupBox_Output.setGeometry(QtCore.QRect(640, 10, 621, 671))
        self.groupBox_Output.setObjectName("groupBox_Output")
        self.Pixmap_Output = QtWidgets.QLabel(self.groupBox_Output)
        self.Pixmap_Output.setGeometry(QtCore.QRect(10, 80, 601, 321))
        self.Pixmap_Output.setText("")
        self.Pixmap_Output.setScaledContents(True)
        self.Pixmap_Output.setObjectName("Pixmap_Output")   
        self.scrollArea_Output = QtWidgets.QScrollArea(self.groupBox_Output)
        self.scrollArea_Output.setGeometry(QtCore.QRect(10, 430, 601, 231))
        self.scrollArea_Output.setWidgetResizable(True)
        self.scrollArea_Output.setObjectName("scrollArea_Output")
        self.scrollAreaWidgetContents_Output = QtWidgets.QWidget()
        self.scrollAreaWidgetContents_Output.setGeometry(QtCore.QRect(0, 0, 599, 229))
        self.scrollAreaWidgetContents_Output.setObjectName("scrollAreaWidgetContents_Output")
        self.plainTextEdit_Output = QtWidgets.QPlainTextEdit(self.scrollAreaWidgetContents_Output)
        self.plainTextEdit_Output.setEnabled(True)
        self.plainTextEdit_Output.setGeometry(QtCore.QRect(0, 0, 601, 241))
        self.plainTextEdit_Output.setReadOnly(True)
        self.plainTextEdit_Output.setPlainText("")
        self.plainTextEdit_Output.setObjectName("plainTextEdit_Output")
        self.scrollArea_Output.setWidget(self.scrollAreaWidgetContents_Output)
        self.horizontalLayoutWidget_Output = QtWidgets.QWidget(self.groupBox_Output)
        self.horizontalLayoutWidget_Output.setGeometry(QtCore.QRect(10, 20, 601, 41))
        self.horizontalLayoutWidget_Output.setObjectName("horizontalLayoutWidget_Output")
        self.horizontalLayout_Output = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_Output)
        self.horizontalLayout_Output.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_Output.setObjectName("horizontalLayout_Output")
        self.lineEdit_Output = QtWidgets.QLineEdit(self.horizontalLayoutWidget_Output)
        self.lineEdit_Output.setObjectName("lineEdit_Output")
        self.horizontalLayout_Output.addWidget(self.lineEdit_Output)
        self.pushButton_Output = QtWidgets.QPushButton(self.horizontalLayoutWidget_Output)
        self.pushButton_Output.setObjectName("pushButton_Output")
        self.horizontalLayout_Output.addWidget(self.pushButton_Output)
        self.comboBox_Output = QtWidgets.QComboBox(self.horizontalLayoutWidget_Output)
        self.comboBox_Output.setObjectName("comboBox_Output")
        self.comboBox_Output.addItems([".txt",".tsv",".csv"])
        self.horizontalLayout_Output.addWidget(self.comboBox_Output)
        self.pushButton_Run = QtWidgets.QPushButton(self.horizontalLayoutWidget_Output)
        self.pushButton_Run.setObjectName("pushButton_Run")
        self.horizontalLayout_Output.addWidget(self.pushButton_Run)
        self.comboBox_Output_Selection = QtWidgets.QComboBox(self.groupBox_Output)
        self.comboBox_Output_Selection.setGeometry(QtCore.QRect(500, 400, 104, 26))
        self.comboBox_Output_Selection.setCurrentText("")
        self.comboBox_Output_Selection.setObjectName("comboBox_Output_Selection")
        self.line = QtWidgets.QFrame(self.groupBox_Output)
        self.line.setGeometry(QtCore.QRect(0, 60, 621, 20))
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")

        self.retranslateUi(Form)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(Form)
        
        self.pushButton_Input.clicked.connect(self.open_file)
        self.pushButton_Recompute.clicked.connect(self.load_sheets)
        self.pushButton_Fingerprints_Preview.clicked.connect(self.get_preview)
        #self.horizontalSlider_Seed.valueChanged.connect(self.scroll_seed)
        self.pushButton_Run.clicked.connect(self.run_script)
        self.pushButton_Output.clicked.connect(self.set_output_path)
        self.comboBox_Output_Selection.currentIndexChanged.connect(self.show_output)
        self.checkbox_fix_column.stateChanged.connect(self.state_changed)
        
        sys.stdout = Stream(newText=self.onUpdateText)
        sys.stderr = Stream(newText=self.onUpdateText)
        
        self.style_reset = self.groupBox_Input.styleSheet()

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.groupBox_Input.setTitle(_translate("Form", "Input Files"))
        self.pushButton_Input.setText(_translate("Form", "Load Data"))
        self.pushButton_Recompute.setText(_translate("Form", "Recompute"))
        self.label_inp_columns.setText(_translate("Form", "Columns to be used for Randomization"))
        self.label_output_columns.setText(_translate("Form", "Further columns to be included in Output"))
        self.label_ref_columns.setText(_translate("Form", "Columns to be used as Label"))
        self.label_Sheets.setText(_translate("Form", "Sheets to be used"))
        
        self.groupBox_Fingerprints_Preview.setTitle(_translate("Form", "Fingerprints"))
        self.pushButton_Fingerprints_Preview.setText(_translate("Form", "Get Fingerprint Preview"))
        
        self.groupBox_Parameters.setTitle(_translate("Form", "Parameters"))
        self.label_Layout.setText(_translate("Form", "Plate Layout"))
        self.label_Blanks.setText(_translate("Form", "Number of Blanks"))
        self.label_Ref1.setText(_translate("Form", "Reference 1"))
        self.label_Ref2.setText(_translate("Form", "Reference 2"))
        
        self.label_num_Ref1.setText(_translate("Form", "#Reference 1"))
        self.label_num_Ref2.setText(_translate("Form", "#Reference 2"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_Fixed), _translate("Form", "Fixed"))
        self.label_min_Ref.setText(_translate("Form", "Minimal #References"))
        self.label_perc_Ref1.setText(_translate("Form", "Percentage Reference 1"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_Relative), _translate("Form", "Relative"))
        
        self.label_Seed.setText(_translate("Form", "Seed"))
        self.label_Settings.setText(_translate("Form", "Settings"))
        self.checkbox_fix_column.setText(_translate("Form", "Group"))
        self.checkbox_optimize.setText(_translate("Form", "Opt."))
        self.checkbox_exclude.setText(_translate("Form", "Excl."))
        
        self.groupBox_Output.setTitle(_translate("Form", "Output"))
        self.pushButton_Output.setText(_translate("Form", "Set Output Path"))
        self.pushButton_Run.setText(_translate("Form", "     Run     "))

    def get_Input(self):
        self.Plate_layout = str(self.comboBox_Layout.currentText())
        if self.Plate_layout == "96":
            self.num_columns = 12
            self.num_rows = 8
        elif self.Plate_layout == "384":
            self.num_columns = 24
            self.num_rows = 16
        self.num_wells = self.num_columns * self.num_rows
        self.Reference_1 = str(self.lineEdit_Ref1.text())
        self.Reference_2 = str(self.lineEdit_Ref2.text())
        self.num_Blanks = int(self.spinBox_Blanks.value())
        self.num_plates = 1
        self.seed = int(self.lineEdit_Seed.text())
        
        print("Input:")
        print("--------------------")
        print("Chosen Plate Layout: {0}".format(self.Plate_layout))
        print("Reference 1: {0}".format(self.Reference_1))
        print("Reference 2: {0}".format(self.Reference_2))
        print("Number of Blanks: {0}".format(self.num_Blanks))
        print("Seed: {0}".format(self.seed))
        
        if self.spinBox_num_Ref1.value()+self.spinBox_num_Ref2.value()!=0:
            self.num_Ref_1 = int(self.spinBox_num_Ref1.value())
            self.num_Ref_2 = int(self.spinBox_num_Ref2.value())
            self.num_wells_to_fill = self.num_wells - self.num_Blanks - self.num_Ref_1 - self.num_Ref_2
            print("Number of Reference 1: {0}".format(self.num_Ref_1))
            print("Number of Reference 2: {0}".format(self.num_Ref_2))
        else:
            self.num_References = int(self.spinBox_min_Ref.value())
            self.percentage_Ref_1 = float(self.doubleSpinBox_perc_Ref1.value())
            print("Minimal Number of References: {0}".format(self.num_References))
            print("Percentage of Reference 1: {0}".format(self.percentage_Ref_1))
            self.num_wells_to_fill = self.num_wells - self.num_Blanks - self.num_References
            
    
    def open_file(self):
        self.file = QtWidgets.QFileDialog(self)
        self.fname, self.mask = self.file.getOpenFileName(self,"Load Data","","Table data (*.csv *.tsv *.txt *.xls *.xlsx *.xlsm)")
        self.lineEdit_Input.setText(self.fname)
        self.in_file = str(self.lineEdit_Input.text())
        self.read_Input()
    
    def set_output_path(self):
        self.outputpath = QtWidgets.QFileDialog(self)
        self.outputpath.setFileMode(QtWidgets.QFileDialog.DirectoryOnly)
        if self.outputpath.exec_() == QtWidgets.QDialog.Accepted:
            self.lineEdit_Output.setText(self.outputpath.selectedFiles()[0])
    
    def get_columns(self):
        self.imp_columns = []
        for button in self.buttongroup_imp_columns.buttons():
            if button.isChecked():
                self.imp_columns.append(button.text())
        self.ref_columns = []
        for button in self.buttongroup_ref_columns.buttons():
            if button.isChecked():
                self.ref_columns.append(button.text())
        self.output_columns = []
        for button in self.buttongroup_output_columns.buttons():
            if button.isChecked():
                self.output_columns.append(button.text())
    
    def get_label_column(self):
        self.label_column = []
        for button in self.ui.buttongroup_cutoff.buttons():
            if button.isChecked():
                self.label_column.append(button.text())
    
    #def scroll_seed(self):
    #    self.lineEdit_Seed.setText(str(self.horizontalSlider_Seed.value()))
    
    def read_Input(self):
        if self.scrollAreaWidgetContents_Sheets.layout():
            QtWidgets.QWidget().setLayout(self.scrollAreaWidgetContents_Sheets.layout()) 
        if self.scrollAreaWidgetContents_imp_columns.layout():
            QtWidgets.QWidget().setLayout(self.scrollAreaWidgetContents_imp_columns.layout())
        if self.scrollAreaWidgetContents_ref_columns.layout():
            QtWidgets.QWidget().setLayout(self.scrollAreaWidgetContents_ref_columns.layout())    
        if self.scrollAreaWidgetContents_output_columns.layout():
            QtWidgets.QWidget().setLayout(self.scrollAreaWidgetContents_output_columns.layout())  
        
        if (self.in_file.endswith(".tsv")) or (self.in_file.endswith(".txt")):
            self.data = pd.read_csv(self.in_file, sep="\t", encoding = "latin")
            self.get_selections_input()
        elif self.in_file.endswith(".csv"):
            self.data = pd.read_csv(self.in_file, sep = ",", encoding = "latin")
            self.get_selections_input()
        elif (self.in_file.endswith(".xlsx")) or (self.in_file.endswith(".xls")) or (self.in_file.endswith(".xlsm")):
            self.data_raw = pd.read_excel(self.in_file, sheet_name=None)
            self.vbox_Sheets = QtWidgets.QVBoxLayout(self)
            self.buttongroup_Sheets = QtWidgets.QButtonGroup(self, exclusive=False)
            for i in list(self.data_raw.keys()):
                self.buttonz = QtWidgets.QCheckBox(str(i),self)
                self.buttonz.toggle()
                self.vbox_Sheets.addWidget(self.buttonz)
                self.buttongroup_Sheets.addButton(self.buttonz)
            self.scrollAreaWidgetContents_Sheets.setLayout(self.vbox_Sheets)
        
    def get_selections_input(self):
        self.vbox_imp_columns = QtWidgets.QVBoxLayout(self)
        self.buttongroup_imp_columns = QtWidgets.QButtonGroup(self, exclusive=False)
        for i in list(self.data.columns):
            self.buttonz = QtWidgets.QCheckBox(str(i),self)
            self.vbox_imp_columns.addWidget(self.buttonz)
            self.buttongroup_imp_columns.addButton(self.buttonz)
        self.scrollAreaWidgetContents_imp_columns.setLayout(self.vbox_imp_columns)
        self.vbox_ref_columns = QtWidgets.QVBoxLayout(self)
        self.buttongroup_ref_columns = QtWidgets.QButtonGroup(self, exclusive=False)
        for i in list(self.data.columns):
            self.buttonz = QtWidgets.QCheckBox(str(i),self)
            self.vbox_ref_columns.addWidget(self.buttonz)
            self.buttongroup_ref_columns.addButton(self.buttonz)
        self.scrollAreaWidgetContents_ref_columns.setLayout(self.vbox_ref_columns)
        self.vbox_output_columns = QtWidgets.QVBoxLayout(self)
        self.buttongroup_output_columns = QtWidgets.QButtonGroup(self, exclusive=False)
        for i in list(self.data.columns):
            self.buttonz = QtWidgets.QCheckBox(str(i),self)
            self.vbox_output_columns.addWidget(self.buttonz)
            self.buttongroup_output_columns.addButton(self.buttonz)
        self.scrollAreaWidgetContents_output_columns.setLayout(self.vbox_output_columns)
        
    def load_sheets(self):
        try:
            list_c = []
            for button in self.buttongroup_Sheets.buttons():
                if button.isChecked():
                    list_c.append(button.text())
            self.data = pd.read_excel(self.in_file, sheet_name=list_c)
            self.data = pd.concat(self.data)
            self.get_selections_input()
        except:   
            self.raiseError("Please select Input File first.")    
    
    def run_script(self):
        self.groupBox_Parameters.setStyleSheet(self.style_reset)
        self.groupBox_Input.setStyleSheet(self.style_reset)
        self.comboBox_Output_Selection.clear()
        try:
            self.get_Input()
            if self.num_wells_to_fill <= 0:
                print("ERROR: Not enough wells to be filled with samples.")
                self.groupBox_Parameters.setObjectName("ColoredGroupBox")
                self.groupBox_Parameters.setStyleSheet("QGroupBox#ColoredGroupBox { border: 1px solid red;}")
                self.raiseError("Parameter selection failed.")
        except:
            print("ERROR: Parameter selection failed.")
            self.groupBox_Parameters.setObjectName("ColoredGroupBox")
            self.groupBox_Parameters.setStyleSheet("QGroupBox#ColoredGroupBox { border: 1px solid red;}")
            self.raiseError("Parameter selection failed.")
        set_seed(self.seed)
        if self.checkbox_fix_column.isChecked():
            try:
                self.get_columns()
                self.get_label_column()
                self.cutoff = int(self.ui.lineEdit_cutoff.text())
                self.data_dict, self.Fingerprints = get_data_dict(self.data, self.label_column, self.imp_columns)
                self.dict_fix_to_label, self.dict_label_to_fix, self.dict_fix_to_fingerprint, self.labels_cut = get_fixed_stats(self.data, self.ref_columns[0], self.label_column[0], self.cutoff)
                self.Fingerprint_IDs, self.Fingerprints_list, self.num_Fingerprint_fixed = get_Fingerprint(self.data_dict, self.Fingerprints, self.imp_columns)
                print("--------------------")
                print("Starting Sample Distribution:")
                self.num_Fingerprint, self.num_plates = get_num_Fingerprint(self.Fingerprints_list, self.num_wells_to_fill, self.dict_fix_to_fingerprint)
                self.Fingerprints_list_full = get_Fingerprints_list_full_fixed(self.Fingerprints_list, self.dict_fix_to_label)
            except:
                print("ERROR: Input selection or Fingerprint generation failed.")
                self.groupBox_Input.setObjectName("ColoredGroupBox")  
                self.groupBox_Input.setStyleSheet("QGroupBox#ColoredGroupBox { border: 1px solid red;}")
                self.raiseError("Input selection or Fingerprint generation failed.")
            if self.spinBox_num_Ref1.value()+self.spinBox_num_Ref2.value()!=0:
                self.Plates_overlap, self.Plates = distribute_samples_fixed(self.Fingerprints_list, self.data_dict, self.num_plates, self.dict_fix_to_label, self.labels_cut, self.num_wells, self.num_Ref_1+self.num_Ref_2, self.num_Blanks) 
                self.Plates_final = distribute_References(self.Plates_overlap, self.Fingerprints_list_full, self.Reference_1, self.Reference_2, self.num_wells, self.num_Blanks, num_Ref_1=self.num_Ref_1, num_Ref_2=self.num_Ref_2, dict_fix_to_label=self.dict_fix_to_label)
            else:
                self.Plates_overlap, self.Plates = distribute_samples_fixed(self.Fingerprints_list, self.data_dict, self.num_plates, self.dict_fix_to_label, self.labels_cut, self.num_wells, self.num_References, self.num_Blanks)             
                self.Plates_final = distribute_References(self.Plates_overlap, self.Fingerprints_list_full, self.Reference_1, self.Reference_2, self.num_wells, self.num_Blanks, percentage_Ref_1=self.percentage_Ref_1, dict_fix_to_label=self.dict_fix_to_label)
            self.data_dict_fixed = get_data_dict_fixed(self.data, self.ref_columns)
            print("--------------------")
            print("Writing Output to: "+str(os.path.join(self.lineEdit_Output.text(),"Out"+self.comboBox_Output.currentText())))
            generate_Output(self.data_dict_fixed, self.Plates_final, self.ref_columns, self.imp_columns, self.num_columns, self.output_columns, self.Reference_1, self.Reference_2, os.path.join(self.lineEdit_Output.text(),"Out"+self.comboBox_Output.currentText()), label_column=self.label_column) 
            self.num_Fingerprint_per_plate = get_Statistics(self.data_dict, self.Fingerprint_IDs, self.Plates, self.imp_columns)
            print_Statistics(self.Fingerprint_IDs, self.Plates_final, self.num_Fingerprint_per_plate, self.Reference_1, self.Reference_2)
            self.data_out, self.Fingerprint_IDs_plate = read_plates(self.imp_columns, os.path.join(self.lineEdit_Output.text(),"Out"+self.comboBox_Output.currentText()))
            plot_PlateLayout(self.data_out, self.ref_columns[0], self.imp_columns, self.Fingerprint_IDs_plate, self.num_rows, self.num_columns, self.num_plates, os.path.join(self.lineEdit_Output.text(),"Plate_"), self.label_column[0]) 
            self.comboBox_Output_Selection.addItems([str(k+1) for k in range(0,self.num_plates)])
            self.comboBox_Output_Selection.setCurrentText("1")
            self.save_Summary_fixed()
            print("--------------------")   
        else:
            try:
                self.get_columns()
                self.data_dict, self.Fingerprints = get_data_dict(self.data, self.ref_columns, self.imp_columns)
                self.Fingerprint_IDs, self.Fingerprints_list, num_Fingerprint = get_Fingerprint(self.data_dict, self.Fingerprints, self.imp_columns)
                self.num_Fingerprint, self.num_plates, self.num_Sample_per_plate, self.num_Samples_overlap = get_Fingerprint_statistics(self.Fingerprints_list, self.num_wells_to_fill, self.num_plates)
                print("--------------------")
                print("Starting Sample Distribution:")
            except:
                print("ERROR: Input selection or Fingerprint generation failed.")
                self.groupBox_Input.setObjectName("ColoredGroupBox")  
                self.groupBox_Input.setStyleSheet("QGroupBox#ColoredGroupBox { border: 1px solid red;}")
                self.raiseError("Input selection or Fingerprint generation failed.")
            self.Plates_overlap = distribute_Samples(self.Fingerprints_list, self.data_dict, self.num_Sample_per_plate, self.num_plates)
            if self.spinBox_num_Ref1.value()+self.spinBox_num_Ref2.value()!=0:
                self.Plates_final = distribute_References(self.Plates_overlap, self.Fingerprints_list, self.Reference_1, self.Reference_2, self.num_wells, self.num_Blanks, num_Ref_1=self.num_Ref_1, num_Ref_2=self.num_Ref_2)
            else:
                self.Plates_final = distribute_References(self.Plates_overlap, self.Fingerprints_list, self.Reference_1, self.Reference_2, self.num_wells, self.num_Blanks, percentage_Ref_1=self.percentage_Ref_1)
            print("--------------------")
            print("Writing Output to: "+str(os.path.join(self.lineEdit_Output.text(),"Out"+self.comboBox_Output.currentText())))
            generate_Output(self.data_dict, self.Plates_final, self.ref_columns, self.imp_columns, self.num_columns, self.output_columns, self.Reference_1, self.Reference_2, os.path.join(self.lineEdit_Output.text(),"Out"+self.comboBox_Output.currentText())) 
            self.num_Fingerprint_per_plate = get_Statistics(self.data_dict, self.Fingerprint_IDs, self.Plates_final, self.imp_columns)
            print_Statistics(self.Fingerprint_IDs, self.Plates_final, self.num_Fingerprint_per_plate, self.Reference_1, self.Reference_2)
            self.data_out, self.Fingerprint_IDs_plate = read_plates(self.imp_columns, os.path.join(self.lineEdit_Output.text(),"Out"+self.comboBox_Output.currentText()))
            plot_PlateLayout(self.data_out, self.ref_columns[0], self.imp_columns, self.Fingerprint_IDs_plate, self.num_rows, self.num_columns, self.num_plates, os.path.join(self.lineEdit_Output.text(),"Plate_"))
            self.comboBox_Output_Selection.addItems([str(k+1) for k in range(0,self.num_plates)])
            self.comboBox_Output_Selection.setCurrentText("1")
            self.save_Summary()
            print("--------------------")
    
    def get_preview(self):
        try:
            self.get_columns()
            if self.checkbox_fix_column.isChecked():
                self.get_label_column()
                data_dict, Fingerprints = get_data_dict(self.data, self.label_column, self.imp_columns)
            else:
                data_dict, Fingerprints = get_data_dict(self.data, self.ref_columns, self.imp_columns)
            Fingerprint_IDs, Fingerprint_list, num_Fingerprint = get_Fingerprint(data_dict, Fingerprints, self.imp_columns)
            for i in range(len(num_Fingerprint)):
                if i==0:
                    self.plainTextEdit_Fingerprints_Preview.setPlainText(str(Fingerprint_IDs[i])+": "+str(num_Fingerprint[i]))
                else:
                    self.plainTextEdit_Fingerprints_Preview.setPlainText(self.plainTextEdit_Fingerprints_Preview.toPlainText()+"\n"+str(Fingerprint_IDs[i])+": "+str(num_Fingerprint[i]))
        except:
            self.raiseError("Please select columns for label and randomization.")
            
    def show_output(self):
        k = self.comboBox_Output_Selection.currentText()
        if self.checkbox_fix_column.isChecked():
            self.pixmap = QtGui.QPixmap(os.path.join(self.lineEdit_Output.text(),"Plate_"+str(k)+"_label.png"))
            self.Pixmap_Output.setPixmap(self.pixmap)            
        else:
            self.pixmap = QtGui.QPixmap(os.path.join(self.lineEdit_Output.text(),"Plate_"+str(k)+".png"))
            self.Pixmap_Output.setPixmap(self.pixmap)
        
    def onUpdateText(self, text):
        cursor = self.plainTextEdit_Output.textCursor()
        cursor.movePosition(QtGui.QTextCursor.End)
        cursor.insertText(text)
        self.plainTextEdit_Output.setTextCursor(cursor)
        self.plainTextEdit_Output.ensureCursorVisible()

    def __del__(self):
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__
        
    def raiseError(self, error):
        msg = QtWidgets.QMessageBox(self)
        msg.setWindowTitle("Error")
        msg.setText("An Error occured")
        msg.setInformativeText(error)
        msg.setIcon(QtWidgets.QMessageBox.Warning)
        msg.setStandardButtons(QtWidgets.QMessageBox.Cancel)
        msg.setDefaultButton(QtWidgets.QMessageBox.Cancel)
        msg.exec_()    
        
    def save_Summary(self):
        with open(os.path.join(self.lineEdit_Output.text(),"Summary.txt"),"w") as file:
            file.write("--------------------\n")
            file.write("Summary\n")
            file.write("--------------------\n")
            file.write("--------------------\n")
            file.write("Input\n")
            file.write("--------------------\n")
            file.write("Chosen Plate Layout: {0}\n".format(self.Plate_layout))
            file.write("Reference 1: {0}\n".format(self.Reference_1))
            file.write("Reference 2: {0}\n".format(self.Reference_2))
            file.write("Number of Blanks: {0}\n".format(self.num_Blanks))
            file.write("Seed: {0}\n".format(self.seed))
            if self.spinBox_num_Ref1.value()+self.spinBox_num_Ref2.value()!=0:
                file.write("Number of Reference 1: {0}\n".format(self.num_Ref_1))
                file.write("Number of Reference 2: {0}\n".format(self.num_Ref_2))
            else:
                file.write("Minimal Number of References: {0}\n".format(self.num_References))
                file.write("Percentage of Reference 1: {0}\n".format(self.percentage_Ref_1))
            file.write("Number of Plates: {0:d}\n".format(self.num_plates))
            file.write("Columns used for Randomisation: ")
            for ind,imp_column in enumerate(self.imp_columns):
                file.write(str(imp_column))
                if ind == len(self.imp_columns)-1:
                    file.write("\n")
                else:
                    file.write(", ")
            file.write("--------------------\n")
            file.write("Fingerprint Statistics\n")
            file.write("--------------------\n")
            for ind,fingerprint in enumerate(self.Fingerprint_IDs):
                file.write(str(fingerprint)+": "+str(self.num_Fingerprint[ind])+"\n")
            file.write("--------------------\n")
            file.write("Plate Statistics\n")
            file.write("--------------------\n")
            for ind_p,plate in enumerate(self.Plates_final):
                file.write("Plate "+str(ind_p+1)+"\n")
                file.write("#"+self.Reference_1+": %d\n" %len(np.where(plate==self.Reference_1)[0]))
                file.write("#"+self.Reference_2+": %d\n" %len(np.where(plate==self.Reference_2)[0]))
                file.write("#Blank: %d\n" %len(np.where(plate=="Blank")[0]))
                for index, fingerprint in enumerate(self.Fingerprint_IDs):
                    file.write("#"+str(fingerprint)+":%d\n" %self.num_Fingerprint_per_plate[index,ind_p])   
                file.write("-----------------------\n")
                
    def save_Summary_fixed(self):
        with open(os.path.join(self.lineEdit_Output.text(),"Summary.txt"),"w") as file:
            file.write("--------------------\n")
            file.write("Summary\n")
            file.write("--------------------\n")
            file.write("--------------------\n")
            file.write("Input\n")
            file.write("--------------------\n")
            file.write("Chosen Plate Layout: {0}\n".format(self.Plate_layout))
            file.write("Reference 1: {0}\n".format(self.Reference_1))
            file.write("Reference 2: {0}\n".format(self.Reference_2))
            file.write("Number of Blanks: {0}\n".format(self.num_Blanks))
            file.write("Seed: {0}\n".format(self.seed))
            if self.spinBox_num_Ref1.value()+self.spinBox_num_Ref2.value()!=0:
                file.write("Number of Reference 1: {0}\n".format(self.num_Ref_1))
                file.write("Number of Reference 2: {0}\n".format(self.num_Ref_2))
            else:
                file.write("Minimal Number of References: {0}\n".format(self.num_References))
                file.write("Percentage of Reference 1: {0}\n".format(self.percentage_Ref_1))
            file.write("Number of Plates: {0:d}\n".format(self.num_plates))
            file.write("Columns used for Randomisation: ")
            for ind,imp_column in enumerate(self.imp_columns):
                file.write(str(imp_column))
                if ind == len(self.imp_columns)-1:
                    file.write("\n")
                else:
                    file.write(", ")
            file.write("--------------------\n")
            file.write("Fingerprint Statistics\n")
            file.write("--------------------\n")
            for ind,fingerprint in enumerate(self.Fingerprint_IDs):
                file.write(str(fingerprint)+": "+str(self.num_Fingerprint_fixed[ind])+"\n")
            file.write("--------------------\n")
            file.write("Plate Statistics\n")
            file.write("--------------------\n")
            for ind_p,plate in enumerate(self.Plates_final):
                file.write("Plate "+str(ind_p+1)+"\n")
                file.write("#"+self.Reference_1+": %d\n" %len(np.where(plate==self.Reference_1)[0]))
                file.write("#"+self.Reference_2+": %d\n" %len(np.where(plate==self.Reference_2)[0]))
                file.write("#Blank: %d\n" %len(np.where(plate=="Blank")[0]))
                for index, fingerprint in enumerate(self.Fingerprint_IDs):
                    file.write("#"+str(fingerprint)+":%d\n" %self.num_Fingerprint_per_plate[index,ind_p])   
                file.write("-----------------------\n")
                
    def state_changed(self):
        if self.checkbox_fix_column.isChecked():
            try:
                self.Form2 = QtWidgets.QWidget()
                self.ui = FixedColumn(self.comboBox_Layout.currentText(), self.data.columns)
                self.ui.setupUi(self.Form2)
                self.Form2.show()
            except:
                self.raiseError("Please select Input File first.")   
            
class FixedColumn(Ui_Form):
    def __init__(self, max_slider, columns):
        super().__init__()
        self.max_slider = int(max_slider)
        self.columns = columns
    def setupUi(self, Form2):
        Form2.setObjectName("Form2")
        Form2.resize(327, 185)
        self.horizontalSlider_cutoff = QtWidgets.QSlider(Form2)
        self.horizontalSlider_cutoff.setGeometry(QtCore.QRect(20, 160, 241, 22))
        self.horizontalSlider_cutoff.setOrientation(QtCore.Qt.Horizontal)
        self.horizontalSlider_cutoff.setObjectName("horizontalSlider_cutoff")
        self.horizontalSlider_cutoff.setMaximum(self.max_slider)
        self.horizontalSlider_cutoff.setSliderPosition(10)
        self.lineEdit_cutoff = QtWidgets.QLineEdit(Form2)
        self.lineEdit_cutoff.setGeometry(QtCore.QRect(270, 160, 41, 21))
        self.lineEdit_cutoff.setObjectName("lineEdit_cutoff")
        self.lineEdit_cutoff.setText("10")
        self.scrollArea_Sheets_cutoff = QtWidgets.QScrollArea(Form2)
        self.scrollArea_Sheets_cutoff.setGeometry(QtCore.QRect(20, 10, 291, 121))
        self.scrollArea_Sheets_cutoff.setWidgetResizable(True)
        self.scrollArea_Sheets_cutoff.setObjectName("scrollArea_Sheets_cutoff")
        self.scrollAreaWidgetContents_cutoff = QtWidgets.QWidget()
        self.scrollAreaWidgetContents_cutoff.setGeometry(QtCore.QRect(0, 0, 289, 119))
        self.scrollAreaWidgetContents_cutoff.setObjectName("scrollAreaWidgetContents_Sheets_cutoff")
        self.scrollArea_Sheets_cutoff.setWidget(self.scrollAreaWidgetContents_cutoff)
        self.label_Sheets_cutoff = QtWidgets.QLabel(Form2)
        self.label_Sheets_cutoff.setGeometry(QtCore.QRect(20, 140, 291, 16))
        self.label_Sheets_cutoff.setObjectName("label_Sheets_cutoff")
        
        self.vbox_cutoff = QtWidgets.QVBoxLayout(self)
        self.buttongroup_cutoff = QtWidgets.QButtonGroup(self, exclusive=True)
        for i in list(self.columns):
            self.buttonz = QtWidgets.QCheckBox(str(i),self)
            self.vbox_cutoff.addWidget(self.buttonz)
            self.buttongroup_cutoff.addButton(self.buttonz)
        self.scrollAreaWidgetContents_cutoff.setLayout(self.vbox_cutoff)
        
        self.horizontalSlider_cutoff.valueChanged.connect(self.scroll_cutoff)

        self.retranslateUi(Form2)
        QtCore.QMetaObject.connectSlotsByName(Form2)

    def retranslateUi(self, Form2):
        _translate = QtCore.QCoreApplication.translate
        Form2.setWindowTitle(_translate("Form2","Fixed Column"))
        self.label_Sheets_cutoff.setText(_translate("Form", "Cutoff for pre-distribution"))
        
    def scroll_cutoff(self):
        self.lineEdit_cutoff.setText(str(self.horizontalSlider_cutoff.value()))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Form = QtWidgets.QWidget()
    ui = Ui_Form()
    ui.setupUi(Form)
    Form.show()
    sys.exit(app.exec_())