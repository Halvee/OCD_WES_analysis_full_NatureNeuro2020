
from collections import OrderedDict
from argparse import ArgumentParser
import pandas

def main():

    # get args
    args = get_args()

    # read xlsx sheet
    df = pandas.read_excel(args.in_xlsx,
                           sheetname=args.in_sheet_name)

    # write to csv
    df.to_csv(path_or_buf=args.out_csv,
              header=True,
              index=False)
       
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
    parser.add_argument('in_xlsx',action='store',
                        type=str, help='name of input xlsx file')
    parser.add_argument('in_sheet_name',action='store',type=str, 
                        help='name of xlsx sheet to subset on')
    parser.add_argument('out_csv',action='store',
                        type=str, help='name of output csv file') 
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    main()
