import numpy as np
import pandas as pd
import hl7
import os
import datetime
import csv
from pathlib import Path
import os.path
import re
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def find_cnv_start_line(search_term):
    """
        Searches for the specified `search_term` in the rows of CSV files with names ending in "Results.csv"
        within the directory specified by the constant QS_INPUT.

        Parameters:
        - search_term (str): The string to search for in the CSV file rows.

        Returns:
        int: The line number (0-indexed) where the first occurrence of `search_term` is found.
             Returns -1 if the term is not found in any of the CSV files.

        result = find_cnv_start_line("desired_term")
        if result != -1:
            print(f"Search term found in CSV file at line {result}.")
        else:
            print("Search term not found in any CSV file.")
        """
    for x in os.listdir(QS_INPUT):
        if x.endswith("Results.csv"):
            a = search_term  # String that you want to search
            num = 0
            with open(x) as f_obj:
                reader = csv.reader(f_obj, delimiter=',')
                for line in reader:
                    if a in str(line):
                        break
                    else:
                        num = num + 1
            return num


def cnv_data_frame(num):
    """
        Reads a CNV (Copy Number Variation) data file, skips the specified number of rows,
        selects relevant columns, and returns a Pandas DataFrame.

        Parameters:
        - num (int): The number of rows to skip in the CNV data file.

        Returns:
        pandas.DataFrame: A DataFrame containing the columns 'Sample_ID', 'Assay ID', 'Call', and 'Gene Symbol'.
        """
    for x in os.listdir(QS_INPUT):
        if x.endswith("Results.csv"):
            cnv_file_name = x
            cnv_df = pd.read_csv(cnv_file_name, skiprows=num, header=0)
            cnv_df = cnv_df[['Sample Name', 'Target', 'CN Predicted']]
            cnv_df.rename(columns={'Sample Name': 'Sample_ID', 'Target': 'Assay ID', 'CN Predicted': 'Call'},
                          inplace=True)
            cnv_df['Gene Symbol'] = 'CYP2D6'
            return cnv_df
        else:
            pass


def find_snp_start_line(search_term):
    """
        Searches for the specified `search_term` in rows of CSV files with names ending in '_Export.csv'
        within the directory specified by the constant QS_INPUT.

        Parameters:
        - search_term (str): The string to search for in the SNP data file rows.

        Returns:
        int: The line number (0-indexed) where the first occurrence of `search_term` is found.
             Returns -1 if the term is not found in any of the SNP data files.
        """
    for x in os.listdir(QS_INPUT):
        if x.endswith('_Export.csv'):
            a = search_term  # String that you want to search
            snp_num = 0
            with open(x, encoding='unicode_escape') as f_obj:
                reader = csv.reader(f_obj, delimiter=',')
                for line in reader:
                    if a in str(line):
                        break
                    else:
                        snp_num = snp_num + 1
            return snp_num


def snp_data_frame(snp_num):
    """
        Reads a SNP (Single Nucleotide Polymorphism) data file, skips the specified number of rows,
        selects relevant columns, and returns a Pandas DataFrame.

        Parameters:
        - snp_num (int): The number of rows to skip in the SNP data file.

        Returns:
        pandas.DataFrame: A DataFrame containing the columns 'Sample_ID', 'Assay ID', 'Call', and 'Gene Symbol'.
        """
    global snp_file_name
    for x in os.listdir(QS_INPUT):
        if x.endswith('_Export.csv'):
            snp_file_name = x
    snp_df = pd.read_csv(snp_file_name, skiprows=snp_num, header=0, encoding='unicode_escape')
    snp_df.dropna(subset=['ROX'], inplace=True)
    snp_df = snp_df[['Sample ID', 'Assay ID', 'Call', 'Gene Symbol']]
    snp_df.rename(columns={'Sample ID': 'Sample_ID'}, inplace=True)
    return snp_df


def snp_cnv_data_frame(cnv_df, snp_df):
    """
        Combines CNV and SNP DataFrames into a single DataFrame.

        Parameters:
        - cnv_df (pandas.DataFrame): DataFrame containing CNV data.
        - snp_df (pandas.DataFrame): DataFrame containing SNP data.

        Returns:
        pandas.DataFrame: A DataFrame combining columns from both CNV and SNP DataFrames.
        """
    full_df = pd.concat([snp_df, cnv_df])
    return full_df


def clean_qs_data(full_df):
    """
    Clean up the raw data from Quantstudios.

    This function takes a DataFrame containing raw data from Quantstudios (QS),
    performs data cleaning operations, and returns the cleaned DataFrame.

    Parameters:
    - full_df (pandas.DataFrame): The input DataFrame containing raw QS data.

    Returns:
    pandas.DataFrame: A cleaned DataFrame with the following operations:
    1. Converts all columns to string type.
    2. Discards rows where 'Sample_ID' contains values specified in 'discard' list.
    3. Updates gene symbols based on predefined replacements.

    Note:
    The 'discard' list and 'gene_symbol_replacements' dictionary are specific to the
    dataset's characteristics and may need adjustment based on the actual data.
    """
    qs_df = full_df
    qs_df = qs_df.astype(str)
    discard = ['CNV Control']
    qs_df = qs_df[~qs_df.Sample_ID.str.contains('|'.join(discard))]

    # update gene symbols from QS output
    gene_symbol_replacements = {
        'CYP2D6,LOC100132273': 'CYP2D6',
        'CYP2D7P1,CYP2D6,LOC100132273': 'CYP2D6',
        'LOC100132273,CYP2D6': 'CYP2D6',
        'NDUFA6-AS1,CYP2D6': 'CYP2D6',
        'CYP2D7P1,CYP2D6': 'CYP2D6',
        'POL3S,VKORC1': 'VKORC1',
        'NUDT15,SUCLA2': 'NUDT15',
        'NHLRC1,TPMT': 'TPMT',
        'CKAP5,F2': 'F2',
        'MTHFR,CLCN6,C1orf167': 'MTHFR'
    }

    # Update gene names using the dictionary
    qs_df['Gene Symbol'] = qs_df['Gene Symbol'].replace(gene_symbol_replacements, regex=True)

    return qs_df


def order_assays(qs_df):
    """
    Ensure assays are in the correct order for gt_pt_translator.

    This function takes a DataFrame with assay data and ensures that the assays are
    ordered based on a custom dictionary. It also performs specific replacements for
    'Gene Symbol' values associated with certain assays.

    Parameters:
    - qs_df (pandas.DataFrame): DataFrame containing assay data with an 'Assay ID' column.

    Returns:
    pandas.DataFrame: A DataFrame with assays ordered based on the custom dictionary.
    'Assay ID' column is renamed to 'Assay_ID', and specific 'Gene Symbol' values are updated.

    Custom Dictionary:
    {'C__25986767_70': 0, 'C__27861809_10': 1, ..., 'C___8726802_20': 46, 'C___1202883_20': 47}

    Note:
    - The custom dictionary defines the order of assays based on their 'Assay ID'.
    - 'Gene Symbol' values for specific assays ('AH1RT7T' and 'AHWR1JA') are updated.
    """

    custom_dict = {'C__25986767_70': 0, 'C__27861809_10': 1, 'C__30634136_10': 2, 'C__27531918_10': 3,
                   'C__30634130_30': 4, 'C____469857_10': 5, 'C__30634128_10': 6, 'Hs00010001_cn': 7,
                   'C___2222771_A0': 8, 'C__11484460_40': 9, 'C__27102414_10': 10, 'C__27102425_10': 11,
                   'C__27102431_D0': 12, 'C__32388575_A0': 13, 'C__32407229_60': 14, 'C__32407232_50': 15,
                    'C__32407243_20': 16, 'C__34816113_20': 17, 'C__34816116_20': 18, 'C_30634117C_K0': 19,
                    'C__26201809_30': 20, 'C__30203950_10': 21, 'C__32287188_10': 22, 'C__30633906_10': 23,
                     'C__25625805_10': 24, 'C__27104892_10': 25, 'C__27859817_40': 26, 'C__32287221_20': 27,
                     'C__25625804_10': 28, 'C__30634132_70': 29, 'C__31983399_10': 30, 'C__16179493_40': 31,
                     'C__30403261_20': 32, 'C_154823200_10': 33, 'C__12091552_30': 34, 'C__30634116_20': 35,
                    'C_____19567_20': 36, 'C__98253221_10': 37, 'AHWR1JA': 38, 'AH1RT7T': 39, 'C___1204093_20': 40,
                    'C___1204091_10': 41, 'C____572770_20': 42, 'C____572771_10': 43, 'C__30026873_10': 44,
                    'C__11975250': 45, 'C___8726802_20': 46, 'C___1202883_20': 47}

    qs_df.sort_values(by=['Assay ID'], key=lambda x: x.map(custom_dict), inplace=True)
    qs_df.rename(columns={'Assay ID': 'Assay_ID'}, inplace=True)
    qs_df.loc[qs_df['Assay_ID'] == 'AH1RT7T', 'Gene Symbol'] = 'APOL1'
    qs_df.loc[qs_df['Assay_ID'] == 'AHWR1JA', 'Gene Symbol'] = 'APOL1'

    return qs_df


def translate_results(qs_df):
    """
    Translate results from Quantstudios (QS) assays to Genotype-Phenotype (GT-PT) information.

    This function utilizes translation tables to map results from QS assays to GT-PT information.
    The QS results are matched with the QS translational table, and the resulting information
    is further merged with the GT-PT translator to provide genotype and phenotype details.

    Parameters:
    - qs_df (pandas.DataFrame): DataFrame containing QS assay results with necessary columns
      like 'Assay_ID', 'Call', 'Sample_ID', 'Gene Symbol', etc.

    Returns:
    pandas.DataFrame: A DataFrame containing translated GT-PT information for each sample and gene.

    Translation Tables:
    - QS_translational_table_df: Translation table for mapping QS results.
    - GT_PT_translator_df: Translation table for mapping GT-PT information.

    Example:
    ```python
    # Assuming 'qs_results' is the DataFrame containing QS assay results
    translated_results = translate_results(qs_results)
    ```

    Note:
    - QS results are matched with the QS translational table based on 'Assay_ID' and 'Call'.
    - Additional processing is performed to handle specific cases ('Hs00010001_cn') and concatenate results.
    - The translated results are then merged with the GT-PT translator based on the 'Call' column.
    - The final DataFrame includes detailed GT-PT information for each sample and gene.
    """

    QS_translational_table_df = pd.read_csv(os.path.join(PGX_Tables, 'QS_Translator.csv'), index_col=False)

    #Geneotype_Phenotype
    GT_PT_translator_df = pd.read_csv(os.path.join(PGX_Tables, 'GT_PT_Translator.csv'),
                                      encoding='unicode_escape', keep_default_na=False)

    # update the results to match QS_translational_table
    QS_results_df = pd.merge(qs_df, QS_translational_table_df, how='left', left_on=['Assay_ID', 'Call'],
                             right_on=['Assay_ID', 'if_Call']).fillna("Not Detected")
    QS_results_df.Then_convert = np.where(QS_results_df.Assay_ID.eq('Hs00010001_cn'), QS_results_df.Call,
                                          QS_results_df.Then_convert)

    # convert QS_results_df to concatenated with probes
    QS_results_df['merge'] = QS_results_df['Assay_ID'] + QS_results_df['Then_convert'].astype(str)
    QS_results_df['Call'] = QS_results_df.groupby(['Sample_ID', 'Gene Symbol'])['merge'].transform(lambda x: ''.join(x))

    # merge the translator with results to give pt gt info
    GT_PT_Results_df = pd.merge(QS_results_df, GT_PT_translator_df, how='left', left_on=['Call'], right_on=['if_'])
    GT_PT_Results_df['if_Call'] = np.where(
        (GT_PT_Results_df['Assay_ID'] == 'Hs00010001_cn') & (GT_PT_Results_df['if_Call'] == 'Not Detected'), "",
        GT_PT_Results_df['if_Call'])

    return GT_PT_Results_df


def full_report_review(GT_PT_Results_df):
    """
    Generate a full review report DataFrame for the reviewer based on Genotype-Phenotype (GT-PT) results.

    This function prepares a review report by extracting relevant columns from the GT-PT results DataFrame.
    Columns such as 'Call', 'if_', and 'merge' are dropped to focus on essential information.
    The resulting DataFrame is then sorted based on 'Sample_ID', 'Gene Symbol', and 'Assay_ID'.
    The column order is also adjusted for better readability in the review report.

    Parameters:
    - GT_PT_Results_df (pandas.DataFrame): DataFrame containing GT-PT results with necessary columns.

    Returns:
    pandas.DataFrame: A review report DataFrame with essential information for each sample and gene.

    Note:
    - The function drops unnecessary columns ('Call', 'if_', 'merge') to focus on essential information.
    - The resulting DataFrame is sorted based on 'Sample_ID', 'Gene Symbol', and 'Assay_ID'.
    - The column order is adjusted to include relevant information such as sample ID, assay ID, probe information, gene symbol,
      translated calls ('if_Call', 'Then_convert'), genotype, phenotype, and activity score.
    """

    GT_PT_Full_Results_df = GT_PT_Results_df.drop(columns=['Call', 'if_', 'merge'])
    GT_PT_Full_Results_df = GT_PT_Full_Results_df.sort_values(by=['Sample_ID', 'Gene Symbol', 'Assay_ID'])

    # Reorder columns
    GT_PT_Full_Results_df = GT_PT_Results_df[['Sample_ID', 'Assay_ID', 'Probe_Information', 'Gene Symbol',
                                              'if_Call', 'Then_convert', 'Genotype', 'Phenotype', 'Activity_Score']]

    return GT_PT_Full_Results_df


def process_results_df(GT_PT_Results_df):
    """
    Process the Genotype-Phenotype (GT-PT) results DataFrame.

    This function drops duplicate rows based on the columns 'Sample_ID', 'Gene Symbol', and 'Genotype',
    sorts the DataFrame by 'Sample_ID' and 'Gene Symbol', and selects specific columns for the processed DataFrame.
    The selected columns include 'Sample_ID', 'Gene Symbol', 'Genotype', 'Phenotype', 'Activity_Score', 'Comment',
    and 'Flag_Note'. The function also renames the 'Gene Symbol' column to 'Gene_Symbol'.

    Parameters:
    - GT_PT_Results_df (pd.DataFrame): The input DataFrame containing Genotype-Phenotype results.

    Returns:
    pd.DataFrame: Processed DataFrame with selected columns and unique rows.
    """

    GT_PT_Results_df.drop_duplicates(subset=['Sample_ID', 'Gene Symbol', 'Genotype'], inplace=True)
    GT_PT_Results_df.sort_values(by=['Sample_ID', 'Gene Symbol'], inplace=True)

    # Create a copy of the DataFrame
    processed_results_df = GT_PT_Results_df.copy()

    # Select the desired columns
    processed_results_df = processed_results_df[['Sample_ID', 'Gene Symbol', 'Genotype', 'Phenotype', 'Activity_Score',
                                                 'Comment', 'Flag_Note']]
    # Rename columns
    processed_results_df.rename(columns={'Gene Symbol': 'Gene_Symbol'}, inplace=True)

    return processed_results_df


def results_for_printing(GT_PT_Results_df):
    """
   Prepare Genotype-Phenotype (GT-PT) results DataFrame for printing.

   This function is a wrapper around the 'process_results_df' function. It returns the processed DataFrame
   obtained from 'process_results_df'.

   Parameters:
   - GT_PT_Results_df (pd.DataFrame): The input DataFrame containing Genotype-Phenotype results.

   Returns:
   pd.DataFrame: Processed DataFrame with selected columns and unique rows.
   """
    return process_results_df(GT_PT_Results_df)


def results_for_HL7(GT_PT_Results_df):
    """
    Prepare Genotype-Phenotype (GT-PT) results DataFrame for HL7 reporting.

    This function is a wrapper around the 'process_results_df' function. It first processes the DataFrame using
    'process_results_df' and then removes rows with NaN values in the 'Genotype' column. The resulting DataFrame
    is suitable for HL7 reporting.

    Parameters:
    - GT_PT_Results_df (pd.DataFrame): The input DataFrame containing Genotype-Phenotype results.

    Returns:
    pd.DataFrame: Processed DataFrame with selected columns, unique rows, and NaN-free 'Genotype' values.
    """
    GT_PT_HL7_Results_df = process_results_df(GT_PT_Results_df)
    GT_PT_HL7_Results_df = GT_PT_HL7_Results_df.dropna(subset=['Genotype'])
    return GT_PT_HL7_Results_df


def alternating_colors(df):
    """
    Generate alternating colors for rows in a DataFrame.

    This function creates a list of alternating colors for rows in a DataFrame, suitable for visualizing
    rows with different colors in a table.

    Parameters:
    - df (pd.DataFrame): The input DataFrame for which alternating colors are generated.

    Returns:
    list: List of alternating colors for rows in the DataFrame.
    """

    alternating_colors = [['white'] * len(df.columns), ['lightgray'] * len(df.columns)] * len(df)
    alternating_colors = alternating_colors[:len(df)]

    return alternating_colors


def print_results_pdf(df1, df2):
    """
    This function creates a PDF file containing tables for reviewing PGX TaqMan Genotyping results.
    It takes two DataFrames (df1 and df2) as input, creates tables for each Sample_ID, and saves the
    resulting PDF with a timestamp.

    Parameters:
    - df1 (pd.DataFrame): The first DataFrame containing PGX TaqMan Genotyping results.
    - df2 (pd.DataFrame): The second DataFrame containing PGX TaqMan Genotyping results.

    Returns:
    None
    """
    dir_path = Path(QS_INPUT)
    sample_files = []
    for path in os.listdir(dir_path):
        if os.path.isfile(os.path.join(path)):
            sample_files.append(path)

    with PdfPages('PGX_TaqMan_Genotyping_Data_Review.pdf') as pdf:
        for Sample_ID in df1['Sample_ID'].unique():
            sample_df1 = df1.loc[df1['Sample_ID'] == Sample_ID]
            sample_df2 = df2.loc[df2['Sample_ID'] == Sample_ID]
            alternating_colors1 = alternating_colors(sample_df1)
            fig, axs = plt.subplots(2)
            # fig.tight_layout()
            axs[0].axis('tight')
            axs[0].axis('off')
            axs[1].axis('tight')
            axs[1].axis('off')
            table1 = axs[0].table(cellText=sample_df1.values,
                                  rowLabels=sample_df1.index,
                                  colLabels=sample_df1.columns,
                                  rowColours=['lightblue'] * len(sample_df1),
                                  colColours=['lightblue'] * len(sample_df1.columns),
                                  cellColours=alternating_colors1,
                                  loc='center')
            table1.auto_set_font_size(False)
            table1.set_fontsize(10)
            table1.auto_set_column_width(col=list(range(len(sample_df1.columns))))
            table1.scale(1, 1.5)

            text = datetime.datetime.now().strftime('%m_%d_%Y_%H_%M_%S')
            axs[0].set_title(
                "\nPGX TaqMan Genotype and Phenotype Result Review\n\n" + text + "\n\n" + str(sample_files) + "\n\n",
                fontweight="bold")

            alternating_colors2 = alternating_colors(sample_df2)
            table2 = axs[1].table(cellText=sample_df2.values,
                                  rowLabels=sample_df2.index,
                                  colLabels=sample_df2.columns,
                                  rowColours=['lightblue'] * len(sample_df2),
                                  colColours=['lightblue'] * len(sample_df2.columns),
                                  cellColours=alternating_colors2,
                                  loc='bottom')
            table2.auto_set_font_size(False)
            table2.set_fontsize(10)
            table2.auto_set_column_width(col=list(range(len(sample_df2.columns))))
            table2.scale(1, 1.5)

            pdf.savefig(fig, bbox_inches='tight')
            plt.close()


def format_for_unity(x):
    """
    Format an integer to ensure it has at least two digits.

    This function takes an integer `x` as input and formats it as a string, ensuring that the resulting string
    has at least two digits. If the input integer is less than 10, a leading zero is added to the string.

    Parameters:
    - x (int): The integer to be formatted.

    Returns:
    str: The formatted string representation of the input integer.
    """
    if (x < 10):
        return "0" + str(x)
    else:
        return str(x)


def get_current_formatted_date():
    """
    Get the current date and time formatted as a string.

    This function retrieves the current date and time using the `datetime` module and formats it as a string.
    The formatted string represents the year, month, day, hour, minute, and second. The formatting ensures that
    each component has at least two digits.

    Returns:
    str: The formatted string representing the current date and time.
    """
    currentDT = datetime.datetime.now()
    if currentDT:
        data = format_for_unity(currentDT.year) + format_for_unity(currentDT.month) + format_for_unity(
            currentDT.day) + format_for_unity(currentDT.hour) + format_for_unity(currentDT.minute) + format_for_unity(
            currentDT.second)

        return data
    return str(currentDT)


def get_ticks(dt):
    return (dt - datetime.datetime(1, 1, 1)).total_seconds() * 10000000


def get_first_obx_index(h):
    idx = 0
    for seg in h:
        if seg[0][0] == 'OBX':
            return idx
        idx += 1
    return -1


def update_comments(comments):
    comments_arr = comments.split("\n")
    i = 1
    h = []
    for comment in comments_arr:
        h.append('NTE|{}|L|{}'.format(i, comment))
        i += 1
    return h


def update_msh_segment(h, current_date):
    if h and h['MSH']:
        for msh_segment in h['MSH']:
            if msh_segment:
                msh_segment[7] = current_date
                msh_segment[8] = ''
                msh_segment[9][0][0] = 'ORU'
                msh_segment[9][0][1] = 'R01'
                msh_segment[10] = get_ticks(datetime.datetime.now())


def update_orc_segment(h):
    if h and h['ORC']:
        for orc_segment in h['ORC']:
            orc_segment[1] = 'RE'


def update_obr_segment(h, current_date):
    if h and h['OBR']:
        for obr_segment in h['OBR']:
            obr_segment[22] = current_date
            obr_segment[25] = 'R'
            obr_segment[27] = '^^^^^R^^'


def update_obx_segment(h):
    if h and h['OBX']:
        for obx_segment in h['OBX']:
            obx_segment[2] = 'ST'
            obx_segment[11] = 'P'
            obx_segment[14] = get_current_formatted_date()
            if len(obx_segment) == 19:
                obx_segment.append(obx_segment[14])
            elif len(obx_segment) >= 19:
                obx_segment[19] = obx_segment[14]


def search_file(file_name_query, path):
    os.chdir(path)
    search_result = glob.glob(file_name_query)
    if not search_result:
        print("Couldn't find a file similar to {}".format(file_name_query))
        return None
    return search_result[0]


# Create LRR and VAR OBX segments
def append_OBX_segments(sample_df, obx_segments):
    """
    Appends OBX (Observation/Result) segments to the given list based on the information in the provided sample DataFrame.

    Parameters:
    - sample_df (DataFrame): A pandas DataFrame containing genetic information for different gene symbols.
    - obx_segments (list): A list to which OBX segments will be appended.

    Returns:
    - list: The updated list of OBX segments.

    This function processes the genetic information in the 'sample_df' DataFrame and creates corresponding HL7 OBX
    segments, which are then appended to the 'obx_segments' list. The function follows specific mappings for gene symbols
    and their associated genotypes, phenotypes, comments, and activity scores. It also generates OBX segments for
    'Gene Studied', 'Genotype Display Name', and various LRR (Laboratory Result Range) segments.

    The function uses the current date as the timestamp for the segments. The 'gene_symbols' list specifies the genes
    that will be considered in the processing.

    Note: The input DataFrame ('sample_df') is expected to have columns like 'Gene_Symbol', 'Genotype', 'Phenotype',
    'Comment', 'Flag_Note', and 'Activity_Score' to extract relevant information.
    """

    obx_segments_len = len(obx_segments)
    i = obx_segments_len + 1
    gi = 'a'
    nte_segments_len = len(obx_segments)
    ni = nte_segments_len + 1
    timestamp = current_date
    gene_symbols = ['CYP2C Cluster', 'CYP2C19', 'CYP2C9', 'CYP2D6', 'CYP3A5', 'CYP4F2', 'SLCO1B1', 'VKORC1']

    # These mappings should match LRR/EAP mappings in epic.
    lrr_gene_obx_mapping = {
        'CYP2C Cluster': {'genotype': '1235884^CYP2C CLUSTER GENOTYPE'},
        'CYP2C19': {'genotype': '1235856^CYP2C19 GENOTYPE', 'phenotype': '123859^CYP2C19 PHENOTYPE', 'comment': 'yes'},
        'CYP2C9': {'genotype': '1235870^CYP2C9 GENOTYPE', 'phenotype': '1235872^CYP2C9 PHENOTYPE',
                   'activity_score': '12358722^CYP2C9 ACTIVITY SCORE', 'comment': 'yes'},
        'CYP2D6': {'genotype': '1235814^CYP2D6 GENOTYPE', 'phenotype': '1238513^CYP2D6 PHENOTYPE',
                   'activity_score': '1235815^CYP2D6 ACTIVITY SCORE', 'comment': 'yes'},
        'CYP3A5': {'genotype': '1235875^CYP3A5 GENOTYPE', 'phenotype': '1235876^CYP3A5 PHENOTYPE', 'comment': 'yes'},
        'CYP4F2': {'genotype': '9865^CYP4F2 GENOTYPE'},
        'SLCO1B1': {'genotype': '1235882^SLCO1B1 GENOTYPE', 'phenotype': '1235883^SLCO1B1 PHENOTYPE'},
        'VKORC1': {'genotype': '1235886^VKORC1 GENOTYPE'},
    }

    # fields needed for OBX segments
    for gene_symbol in gene_symbols:
        if gene_symbol in sample_df['Gene_Symbol'].values:
            gene_name = sample_df.loc[sample_df['Gene_Symbol'] == gene_symbol, 'Gene_Symbol'].values[0]
            genotype = sample_df.loc[sample_df['Gene_Symbol'] == gene_symbol, 'Genotype'].values[0]
            phenotype = sample_df.loc[sample_df['Gene_Symbol'] == gene_symbol, 'Phenotype'].values[0]
            comment = sample_df.loc[sample_df['Gene_Symbol'] == gene_symbol, 'Comment'].values[0]
            flag = sample_df.loc[sample_df['Gene_Symbol'] == gene_symbol, 'Flag_Note'].values[0]
            activity_score = sample_df.loc[sample_df['Gene_Symbol'] == gene_symbol, 'Activity_Score'].values[0]

            # VAR required segments (not currently reporting phenotype in var)
            obx_segments.append('OBX|{}|CWE|48018-6^Gene Studied^LN|4{}|^{}^HGNC'.format(i, gi, gene_name))
            i += 1
            obx_segments.append('OBX|{}|ST|12303110046^Genotype Display Name^LN|4{}|{}'.format(i, gi, genotype))
            i += 1

            # VAR for activity score if required -
            if 'activity_score' in lrr_gene_obx_mapping[gene_symbol] and activity_score is not None:
                if '^' in activity_score:
                    obx_segments.append('OBX|{}|NR|123031010474^ACTIVITY SCORE^LN|4{}|{}'.format(i, gi, activity_score))
                else:
                    obx_segments.append('OBX|{}|NM|123031010474^ACTIVITY SCORE^LN|4{}|{}'.format(i, gi, activity_score))
                i += 1
            gi = chr(ord(gi) + 1)

            #LRR segments:
            # LRR for genotype
            obx_segments.append('OBX|{}|ST|{}|123050000|{}||-||||P|||{}|||2|UFHPL GatorSeq|{}'.format(
                i, lrr_gene_obx_mapping[gene_symbol].get('genotype', ''), genotype, timestamp, timestamp))
            i += 1

            # LRR for phenotype if required
            if 'phenotype' in lrr_gene_obx_mapping[gene_symbol] and phenotype is not None:
                obx_segments.append('OBX|{}|ST|{}|123050000|{}||-|{}|||P|||{}|||2|UFHPL GatorSeq|{}'.format(
                    i, lrr_gene_obx_mapping[gene_symbol].get('phenotype', ''), phenotype, flag, timestamp, timestamp))
                i += 1

            # LRR for Comment if required
            if 'comment' in lrr_gene_obx_mapping[gene_symbol] and comment is not None:
                obx_segments.append('NTE|{}|L|{}'.format(ni, comment))
                ni += 1

            # LRR for activity score if required
            if 'activity_score' in lrr_gene_obx_mapping[gene_symbol] and activity_score is not None:
                obx_segments.append('OBX|{}|ST|{}|123050000|{}||-|{}|||P|||{}|||2|UFHPL GatorSeq|{}'.format(
                    i, lrr_gene_obx_mapping[gene_symbol].get('activity_score', ''), activity_score, flag, timestamp,
                    timestamp))
                i += 1

    return obx_segments


def write_hl7(outfile):
    with open(outfile, 'w', encoding='utf-8') as f:
        f.write(str(h))
        f.write("\n")
        f.write(obx_segments_string)


if __name__ == "__main__":
    script_dir = os.path.dirname(__file__)
    PGX_Tables = os.path.join(script_dir, 'PGX_Tables')
    HL7_Results = os.path.join(script_dir, 'HL7_Results')
    HL7_Orders = os.path.join(script_dir, 'HL7_Orders')
    QS_INPUT = os.path.join(script_dir, 'QS_Data')
    os.chdir(QS_INPUT)
    
    num = find_cnv_start_line("Sample Name")
    cnv_df = cnv_data_frame(num)
    snp_num = find_snp_start_line("Assay Name")
    snp_df = snp_data_frame(snp_num)
    full_df = snp_cnv_data_frame(cnv_df, snp_df)
    qs_df = clean_qs_data(full_df)
    qs_df = order_assays(qs_df)
    GT_PT_Results_df = translate_results(qs_df)

    # Full data review report
    GT_PT_Full_Results_df = full_report_review(GT_PT_Results_df)
    GT_PT_Full_Results_df = GT_PT_Full_Results_df.replace(np.nan, '', regex=True)
    GT_PT_HL7_Results_df = results_for_HL7(GT_PT_Results_df)
    unique_samples_df = GT_PT_HL7_Results_df.drop_duplicates(subset=['Sample_ID'])

    # sample results summary review report
    GT_PT_printing_Results_df = results_for_printing(GT_PT_Results_df)
    GT_PT_Review_df = GT_PT_printing_Results_df[['Sample_ID', 'Gene_Symbol', 'Genotype', 'Phenotype', 'Activity_Score']]
    print_results_pdf(GT_PT_Review_df, GT_PT_Full_Results_df)

    #retreive container IDs from HL7 orders and append data for matching sample.
    hl7_folder_path = HL7_Orders

    # Convert the HL7 message to a list of lines
    def extract_sample_ids(hl7_message_order):
        hl7_lines = hl7_message_order.split('\n')
        return [re.search(r'OBR\|1\|\d+\^\w+\|(\d+)\^', line).group(1) for line in hl7_lines if 'OBR|1' in line]


    for filename in os.listdir(hl7_folder_path):
        if filename.endswith(".txt"):  # Adjust the file extension as needed
            file_path = os.path.join(hl7_folder_path, filename)
            with open(file_path, 'r') as file:
                hl7_message_order = file.read()
                extracted_ids = extract_sample_ids(hl7_message_order)

                for sample_id in extracted_ids:
                    # Check if Sample_ID exists in GT_PT_Results_df
                    if sample_id in GT_PT_Results_df['Sample_ID'].values:
                        sample_df = GT_PT_HL7_Results_df.loc[GT_PT_HL7_Results_df['Sample_ID'] == sample_id]
                        msg_unix_fmt = hl7_message_order.replace("\n", "\r")
                        h = hl7.parse(msg_unix_fmt)
                        obx_segments = []
                        current_date = get_current_formatted_date()
                        append_OBX_segments(sample_df, obx_segments)
                        update_msh_segment(h, current_date)
                        update_orc_segment(h)
                        update_obr_segment(h, current_date)
                        obx_segments_string = "\n".join(obx_segments)

                        # hl7 message paths
                        out_file_path = HL7_Results + '/hl7-GatorPGX-{}-output.txt'.format(sample_id)
                        write_hl7(out_file_path)