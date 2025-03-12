#!/usr/bin/env Python3.10
# GWAS_QC_Code.py

"""
Abigail Parakoyi
July 18th, 2024
:
Build a script to standardize and quality check GWAS data.
"""
import argparse
import yaml
import gwaslab as gl


def main():
    """ main function """

    # create cores variable to use all throughout
    cores_available = '20'
    num_cores = int(cores_available)

    # take command line arguments
    args = get_cli_args()

    # Take in path to GWAS and config file
    infile = args.infile
    config = args.config

    # should take in path to output file and not be hard coded
    outfile = args.outfile

    # opening config file
    with open(config, 'r') as f:
        config_output = yaml.safe_load(f)

    # defining attributes in config file
    col_snpid = config_output['col_names']['snpid']
    col_rsid = config_output['col_names']["rsid"]
    col_chrom = config_output['col_names']['chromosome']
    col_bp = config_output['col_names']['base_pair_position']
    col_ea = config_output['col_names']['effect_allele']
    col_nea = config_output['col_names']['non_effect_allele']
    col_eaf = config_output['col_names']['effect_allele_frequency']
    col_maf = config_output['col_names']['minor_allele_frequency']
    col_samp = config_output['col_names']['sample_size']
    col_beta = config_output['col_names']['beta']
    col_se = config_output['col_names']['standard_error']
    col_p = config_output['col_names']['p_value']
    col_or = config_output['col_names']['odds_ratio']
    col_i2 = config_output['col_names']['I2']
    col_cases = config_output['col_names']['N_CASE']
    col_control = config_output['col_names']['N_CONTROL']
    anc = config_output['ancestry']

    # defining ref_sequence and ancestry path
    ref_seq_path = "/PHShome/ap1231/lsf/GWAS_QC/ref_alleles_files/hg38.fa"
    anc_paths = {
        'EUR': "/PHShome/ap1231/lsf/GWAS_QC/ref_alleles_files/1kg_eur_hg38",
        'EAS': "/PHShome/ap1231/lsf/GWAS_QC/ref_alleles_files/1kg_eas_hg38",
        'SAS': "/PHShome/ap1231/lsf/GWAS_QC/ref_alleles_files/1kg_sas_hg38",
        'AMR': "/PHShome/ap1231/lsf/GWAS_QC/ref_alleles_files/1kg_amr_hg38",
        'AFR': "/PHShome/ap1231/lsf/GWAS_QC/ref_alleles_files/1kg_afr_hg38",
        'PAN': "/PHShome/ap1231/lsf/GWAS_QC/ref_alleles_files/1kg_pan_hg38"
    }

    # load in the GWAS file
    ofile = open_GWAS_file(file=infile, col_snpid=col_snpid, col_rsid=col_rsid, col_chrom=col_chrom, col_bp=col_bp,
                           col_ea=col_ea, col_nea=col_nea, col_eaf=col_eaf, col_maf=col_maf, col_samp=col_samp,
                           col_beta=col_beta, col_se=col_se, col_p=col_p, col_or=col_or, col_i2=col_i2,
                           col_control=col_control, col_cases=col_cases)

    # plot qq plot pre QC data
    create_QQ_plot(file=ofile, name="pre_QC_QQ_plot")

    # run basic check standardization on data
    cdata = check_basics(checked_file=ofile)
    gl.dump_pickle(cdata, outfile, overwrite=True)

    # graph qq plot post QC
    create_QQ_plot(cdata, name="post_QC_QQ_plot")

    # lift over of GWAS file to build 38
    lfile = change_build(file=cdata, cores=num_cores)
    gl.dump_pickle(lfile, outfile, overwrite=True)

    # check and or convert beta column
    bfile = convert_or_to_beta(file=cdata)
    gl.dump_pickle(bfile, outfile, overwrite=True)

    # Check and or convert p column
    pfile = convert_mlogp_to_pval(file=bfile)
    gl.dump_pickle(pfile, outfile, overwrite=True)

    # check and or create maf column
    mafile = give_maf(file=pfile)
    gl.dump_pickle(mafile, outfile, overwrite=True)

    # get path to vcf and ref file
    anc_path = get_anc_path(anc=anc, paths=anc_paths)

    # harmonize data
    hfile = harmonize_data(file=mafile, anc=anc_path, cores=num_cores, seq_path=ref_seq_path)
    gl.dump_pickle(hfile, outfile, overwrite=True)

    # calculate and graph gene inflation
    #gifile = calc_gen_inflation(file=hfile, anc=anc_path, cores=num_cores)
    #graph_gen_inflation(fileplot=gifile)

    # save log and Standardized GWAS file
    hfile.data.to_csv(outfile, sep='\t', index=False)


# open and load in file according to yaml config file
def open_GWAS_file(file, col_snpid, col_rsid, col_chrom, col_bp, col_ea, col_nea, col_eaf,
                   col_maf, col_samp, col_beta, col_se, col_p, col_or, col_i2, col_cases, col_control):
    mysumstats = gl.Sumstats(file, snpid=col_snpid, rsid=col_rsid, chrom=col_chrom, pos=col_bp, ea=col_ea, nea=col_nea,
                             eaf=col_eaf, maf=col_maf, n=col_samp, beta=col_beta, se=col_se, p=col_p, OR=col_or,
                             i2= col_i2, ncase=col_cases, ncontrol=col_control)
    return mysumstats


# Create a QQ plot to show the distribution of the data pre- and post-QC
# checking normalization of data and adding to log
def create_QQ_plot(file, name):
    graph = file.plot_mqq(skip=3, cut=20, mode="qq", anno=True, check=False, save=name, verbose=False,
                          stratified=True, save_args={"dpi": 200, "facecolor": "white"})
    return graph


# Run the file through a basic check (sanity, syntax) of the data set
# then update onto the output file
def check_basics(checked_file):
    # gl.load_pickle(checked_file)
    checked_file.basic_check(fixid_args={'fixchrpos': True, 'fixsep': True},
                             removedup_args={'mode': "md", 'keep': "first", 'remove': True})
    return checked_file


# convert the build of the passed GWAS file to build 38
# double check this one
def change_build(file, cores):
    #bd_file = file.infer_build()
    #if bd_file "hg38"
    file.liftover(n_cores=int(cores), from_build="19", to_build="38", remove=False)
    return file


# Check and convert or create beta column within data set
# should return either full data set that includes a beta column or an error message.
def convert_or_to_beta(file):
    if "BETA" in file.data:
        print("BETA value column already present in file")
    elif "OR" in file.data:
        file.fill_data(to_fill=["BETA"])
    else:
        print("Error! OR column is required to calculate BETA and dataset does not have an OR Column.")
    return file


# Check and convert p-values within the data set
# should return updated or the same p-values column in full datat set in correct format
def convert_mlogp_to_pval(file):
    if "P" in file.data:
        print("P value column already present in file")
    elif "MLOG10P" in file.data:
        file.fill_data(to_fill=["P"])
    else:
        print("Error! MLOG10P column is required to calculate P and dataset does not have an MLOG10P Column.")
    return file


# Check and create minor allele frequency column in data set
# should return updated or the same maf column in full datat set in correct format
def give_maf(file):
    if "MAF" in file.data:
        print("MAF column already present in file")
    elif "EAF" in file.data:
        file.fill_data(to_fill=["MAF"])
    else:
        print("Error! EAF column is required to calculate MAF and dataset does not have an EAF Column.")
    return file


# corrected function to get ancestry file path 
def get_anc_path(anc, paths):
    if anc in paths:
        return paths[anc]
    else:
        return "File path not found the given ancestry"


def harmonize_data(file, anc, cores, seq_path):
    file.harmonize(basic_check=False,
                   n_cores=int(cores),
                   ref_seq=seq_path,
                   ref_infer=anc,
                   ref_alt_freq="AF")
    return file


# Calculate the genomic inflation
# Should return graph of the genomic inflation and add to the log
def calc_gen_inflation(file, anc, cores):
    file.check_af(ref_infer=anc,
                  ref_alt_freq="AF",
                  n_cores=int(cores))
    return file


# Graph genomic inflation
def graph_gen_inflation(fileplot):
    fileplot.plot_daf(threshold=0.12)
    return fileplot


# save the log files and give as output
# should return the completely ran log file
def save_log(file):
    file.log.save()
    return file


def get_cli_args():
    """
    Access the code through the command line
    :return: Instance of argparse argument
    """
    parser = argparse.ArgumentParser(description='Provide path to  GWAS yaml config files  to '
                                                 'QC and standardize')
    parser.add_argument('-i', '--infile',
                        dest='infile',
                        type=str,
                        help='Path to GWAS files to open',
                        required=True)
    parser.add_argument('-c', '--config',
                        dest='config',
                        type=str,
                        help='Path to yaml config files to open',
                        required=False)
    parser.add_argument('-o', '--outfile',
                        dest='outfile',
                        type=str,
                        help='Path to txt file to write',
                        required=True)
    return parser.parse_args()


if __name__ == '__main__':
    main()


