import pandas as pd

class PandasHandler:

    def __init__(self):
        print("This is Constructor")
        
    def __del__(self):
        print("This is Destructor")

    @staticmethod
    def getSheetNumber_Excel(filePathExcel=None):
        bk = pd.ExcelFile(filePathExcel)
        return len(bk.sheet_names)

    @staticmethod
    def getSheetNames_Excel(filePathExcel=None):
        bk = pd.ExcelFile(filePathExcel)
        return bk.sheet_names

    @staticmethod
    def readAllSheets_Excel(filePathExcel=None):
        numSheets = PandasHandler.getSheetNumber_Excel(filePathExcel)
        dfSet = []
        for i in range(numSheets):
            df  = pd.read_excel(filePathExcel, sheet_name=i, index_col=None)
            dfSet.append(df)
        return dfSet
