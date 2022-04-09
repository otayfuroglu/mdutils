import MDAnalysis
import numpy as np
import os
import pandas as pd
import itertools
import argparse


def getDataFromXVG(file_name):
    """
    args
        file name: xvg file
    returns
        2D numpy array
    """
    aux = MDAnalysis.auxiliary.XVG.XVGReader("%s/%s" %(xvgdir, file_name))
    data = np.vstack([row.data for row in aux])
    return data


def padNan(data_ls):
    # pedding nan 
    max_len = max([arr.shape for arr in data_ls])
    data_pad_ls = []
    for data in data_ls:
        arr = np.zeros(max_len)
        arr[:] = np.nan
        for i, row in enumerate(data):
            arr[i] = row
        data_pad_ls.append(arr)
    return data_pad_ls


def getAvgStep2(groupedXvgFile):

    data_ls = [getDataFromXVG(xvgFile) for xvgFile in groupedXvgFile]
    data_pad_ls = padNan(data_ls)
    data_mean = np.nanmean(np.array(data_pad_ls), axis=0)

    df = pd.DataFrame()
    df["Time"] = pd.Series(data_mean[:, 0])
    df["Avg"] = pd.Series(data_mean[:, 1])

    if not os.path.exists(workdir):
        os.mkdir(workdir)
    csv_name = groupedXvgFile[0][:-5]
    df.to_csv("%s/%savg.csv" %(workdir, csv_name), index = False)

    return data_mean[:, 1]


def autoDataFinder():

    files = os.listdir(xvgdir)
    xvgFiles = [fl for fl in files if ".xvg" in fl and fl[-5:-4].isnumeric() and not "cluster" in fl]
    xvgFiles = sorted(xvgFiles) # sorted for grouped func works correctly.

    df = pd.DataFrame()
    groupedXvgFiles = [list(i) for j, i in itertools.groupby(xvgFiles, lambda a: a.split('_')[0])]

    data_ls = [getDataFromXVG(xvgFile) for xvgFile in groupedXvgFiles[0]]

    data_pad_ls = padNan(data_ls)
    data_mean = np.nanmean(np.array(data_pad_ls), axis=0)
    df["Time"] = pd.Series(data_mean[:, 0])

    for i, groupedXvgFile in enumerate(groupedXvgFiles):
        print("Averege for: ", groupedXvgFile)
        data = getAvgStep2(groupedXvgFile)

        df.loc[:, groupedXvgFile[0][:-6]] = pd.Series(data)
    df.to_csv("%s/all_avg.csv" %workdir, index = False)

if __name__ == "__main__":
    workdir = "avgs"

    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-auto", type=str, required=True, help="")
    parser.add_argument("-xvg_dir", type=str, required=True, help="")
    parser.add_argument("-groupedXvgBase", "--groupedXvgBase", type=str, required=False, help="2D file name list")
    args = parser.parse_args()

    xvgdir = args.xvg_dir
    groupedXvgBase = args.groupedXvgBase
    print(groupedXvgBase)

    if args.auto == "on":
        autoDataFinder()
    else:
        files = os.listdir(xvgdir)
        groupedXvgFile = [fl for fl in files if ".xvg" in fl and groupedXvgBase in fl]

        print("Averege for: ", groupedXvgFile)
        data = getAvgStep2(groupedXvgFile)

        if os.path.exists("%s/all_avg.csv" %workdir):
            df = pd.read_csv("%s/all_avg.csv" %workdir, index_col = False)
            df.loc[:, groupedXvgFile[0][:-6]] = pd.Series(data[:, 1])
            df.to_csv("%s/all_avg.csv" %workdir, index = False)
        else:
            print("all_avg.csv file not FOUND!")
