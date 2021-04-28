
from collections import OrderedDict
from argparse import ArgumentParser
import pandas

def main():

    # get args
    args = get_args()

    # read legend
    legend, legend_full = read_legend_txt(args.legend_txt)

    # init an ExcelWriter Obj for writing a single xlsx file
    ew = pandas.ExcelWriter(args.out_xlsx, engine='xlsxwriter')

    # get workbook obj from ExcelWriter obj
    workbook  = ew.book

    # if legend txt defined, write first sheet as a 'Legend',
    # fill in this sheet at the very end of script.
    if args.legend_txt != None:
        cell_format = workbook.add_format()
        cell_format.set_align('left')
        cell_format.set_text_wrap()
        bold = workbook.add_format({'bold': True})
        left = workbook.add_format({'align': 'left', 'text_wrap': True})
        tablenames = list(legend.keys())
        legend_sheet = workbook.add_worksheet("Legend")
        legend_sheet.set_column(0, 0, 100)

    # setup workbook obj for center-aligned text
    format = workbook.add_format()
    format.set_align('center')
    format.set_align('vcenter')

    # init dictionary for storing dfs
    dfs = dict()

    for i in range(len(args.sheet_csv_names)):

        # get sheet name, csv name
        sheet_csv_name = args.sheet_csv_names[i].split("%")
        sheet_name = sheet_csv_name[0]
        csv_name = sheet_csv_name[1]

        # read csv into dataframe
        df = pandas.read_csv(csv_name, sep=args.input_table_delim)
        print("read " + str(df.shape[0]) + " rows from " + \
              csv_name + ", placing into sheet '" + sheet_name + "'.")

        # store df to dfs dictionary object
        dfs[sheet_name] = df

        # if cols_namechange defined, change column names
        cols_namechange={}
        if args.cols_namechange != None:
            # remove '_variant' from Effect column
            effs = list(df["Effect"])
            for j in range(len(effs)):
                effs[j] = effs[j].replace("_variant","")
            df["Effect"] = effs
            for oldnew in args.cols_namechange.split(","):
                [old, new] = oldnew.split(":")
                cols_namechange[old] = new
            columns_new = []
            for j in range(len(df.columns)):
                if df.columns[j] in cols_namechange:
                    columns_new.append(cols_namechange[ df.columns[j] ])
                else:
                    columns_new.append(df.columns[j])
            df.columns = columns_new

        # write to ExcelWriter object
        df.to_excel(ew, 
                    sheet_name=sheet_name, 
                    startrow=1,
                    index=False)

    # auto-set column widths in each sheet.
    # adapted from :
    # https://stackoverflow.com/questions/17326973/is-there-a-way-to-auto-adjust-excel-column-widths-with-pandas-excelwriter
    center = workbook.add_format({'align': 'center'})
    for sheet_name in dfs.keys():
       df = dfs[sheet_name]
       for idx, col in enumerate(df.columns):  # loop through all columns
           series = df[col]
           max_len = max((
           series.astype(str).map(len).max(),  # len of largest item
           len(str(series.name))  # len of column name/header
           )) + 1  # adding a little extra space
           # center-align the text using previously created workbook obj
           ew.sheets[sheet_name].set_column(idx, idx, center)
           # set column width
           ew.sheets[sheet_name].set_column(idx, idx, max_len)  

    # get workbook obj from ExcelWriter obj, for each sheet..
    workbook  = ew.book
    i = 0
    while i < len(workbook.formats):
        # center-align text in all cells
        workbook.formats[i].set_align('center')
        i += 1

    # for each sheet ..
    cell_format = workbook.add_format()
    cell_format.set_align('left')
    bold = workbook.add_format({'bold': True})
    for sheet_name in dfs.keys():    
        # write empty string to first cell that to left-align it
        if args.legend_txt == None:
            ew.sheets[sheet_name].write_blank(0, 0, None, cell_format)
        else:
            ew.sheets[sheet_name].write_rich_string(0, 0, 
                                                    bold, sheet_name+". ",
                                                    legend[sheet_name], cell_format)      

    # if legend provided, fill in the Legend sheet as final step
    if args.legend_txt != None:
        bold = workbook.add_format({'bold': True})
        left = workbook.add_format({'align': 'left', 'text_wrap': True})
        for i in range(len(tablenames)):
            legend_sheet.write_rich_string(i, 0,
                                           bold, tablenames[i]+". ",
                                           legend[tablenames[i]] + " " + \
                                           legend_full[tablenames[i]], left)

    # close ExcelWriter connection
    ew.save()
        
    return

def read_legend_txt(legend_txt_file):
    tables_legend=OrderedDict()
    tables_legend_full=OrderedDict()
    if legend_txt_file == None:
        return tables_legend, tables_legend_full
    fh = open(legend_txt_file,"r")
    current_table_name=None
    for line in fh:
        line = line.rstrip()
        if line[:2] == "##":
            info = line[2:].split(":")
            table_name = info[0]
            table_desc = info[1]
            tables_legend[table_name] = table_desc
            tables_legend_full[table_name] = ""
            current_table_name=table_name
        else:
            line=line+" "
            tables_legend_full[current_table_name]+=line
    fh.close()
    return tables_legend, tables_legend_full

def get_args():
    parser = ArgumentParser(description='')
    parser.add_argument('--cols-keep',action='store',default=None,type=str,
                        help='name of columns to subset on and keep.')
    parser.add_argument('--cols-namechange',action='store',default=None,type=str,
                        help='name changes to apply to cols (ex: ' + \
                        'old1:new1,old2:new2)')
    parser.add_argument('--ids-replace-csv',action='store',default=None,type=str,
                        help='csv with old->new sample IDs (col1->col2)')
    parser.add_argument('--input-table-delim',action='store',default=',',type=str,
                        help='delimiter char for input tables')
    parser.add_argument('--legend-txt', action='store', default=None,
                        type=str, help='name of legend text file')
    parser.add_argument('--out-xlsx',action='store',default='merged_csvs.xlsx',
                        type=str, help='name of output xlsx file')
    parser.add_argument('sheet_csv_names',action='store',type=str, nargs='+')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    main()
