import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser

os.environ.setdefault("RPY2_CFFI_MODE", "ABI")

import rpy2.robjects.packages as rpackages
import rpy2.robjects as ro
import pickle



def zero_one_to_bool(s: str) -> bool:
    try:
        v = int(s)
    except ValueError:
        raise ArgumentTypeError("must be 0 or 1")
    if v not in (0, 1):
        raise ArgumentTypeError("must be 0 or 1")
    return bool(v)

def none_if_empty(x):
    return None if isinstance(x, str) and x == "" else x

# the R object is composed of nested lists

def convert_matrix_to_dataframe(matrix):
    array = np.array(matrix)
    rownames = ro.r['rownames'](matrix)
    colnames = ro.r['colnames'](matrix)
    
    if rownames.rclass[0] == 'NULL':
        rownames = None
    if colnames.rclass[0] == 'NULL':
        colnames = None
        
    df = pd.DataFrame(array, index=rownames, columns=colnames)
    return df

def convert_named_vector_to_series(named_vector):
    names = ro.r['names'](named_vector)
    if names.rclass[0] != 'NULL':
        return pd.Series(named_vector, index=names)
    else:
        return list(named_vector)

def parse_element(element):
    if isinstance(element, ro.vectors.ListVector):
        return {k: parse_element(element.rx2(k)) for k in element.names}
    elif isinstance(element, ro.vectors.Matrix):
        return convert_matrix_to_dataframe(element)
    elif isinstance(element, ro.vectors.Vector):
        return convert_named_vector_to_series(element)
    else:
        return element

def parse_r_object(r_object):
    return {k: parse_element(r_object.rx2(k)) for k in r_object.names}

def save_data_to_pickle(file_path, data):
    with open(file_path, 'wb') as file:
        pickle.dump(data, file)

## run MixTCRviz function from the MixTCRviz package


if __name__ == '__main__':
    parser = ArgumentParser()

    # Core I/O (strings -> pass None to become NULL in R)
    parser.add_argument('-i', '--input1', default=None)
    parser.add_argument('-o', '--output_path', default=None)      # maps to output.path
    parser.add_argument('--input2', default="")                   # -> NULL if ""
    parser.add_argument('--baseline_file', default="")            # -> NULL if ""

    # Chain & protocol
    parser.add_argument('--chain_list_output', default="AB")      # chain
    parser.add_argument('--seq_protocol', default="Default")      # seq.protocol

    # Interactivity / plotting style
    parser.add_argument('--interactive_plots', type=zero_one_to_bool, default=False)  # interactive.plots
    parser.add_argument('--plot', type=zero_one_to_bool, default=True)
    parser.add_argument('--plot_cdr12_motif', type=zero_one_to_bool, default=False)   # plot.cdr12.motif
    parser.add_argument('--plot_oneline', type=int, default=0)                         # integer mode
    parser.add_argument('--plot_logo_length', type=zero_one_to_bool, default=False)    # plot.all.length
    parser.add_argument('--plot_cdr3_norm', type=int, default=0)                       # integer mode
    parser.add_argument('--plot_VJ_switch', type=int, default=1)                       # plot.VJ.switch
    parser.add_argument('--plot_modelsCombined', type=zero_one_to_bool, default=False) # plot.modelsCombined
    parser.add_argument('--plot_sd', type=zero_one_to_bool, default=True)              # plot.sd
    parser.add_argument('--plot_title', type=zero_one_to_bool, default=True)           # plot.title
    parser.add_argument('--set_title', default="")                                     # set.title -> NULL if ""

    # Labels / thresholds
    parser.add_argument('--label_neg', type=zero_one_to_bool, default=False)           # label.neg
    parser.add_argument('--label_diag', type=float, default=0.3)                       # label.diag
    parser.add_argument('--label_min_fr_input1', type=float, default=0.05)             # label.min.fr.input1
    parser.add_argument('--label_min_fr_input2', type=float, default=0.05)             # label.min.fr.input2
    parser.add_argument('--ZscoreVJ_thresh', type=float, default=0.0)                  # ZscoreVJ.thresh
    parser.add_argument('--FoldChangeVJ_thresh', type=float, default=1.25)             # FoldChangeVJ.thresh

    # Gene/allele handling
    parser.add_argument('--use_allele', type=zero_one_to_bool, default=False)          # use.allele
    parser.add_argument('--correct_gene_names', type=zero_one_to_bool, default=True)   # correct.gene.names
    parser.add_argument('--use_mouse_strain', type=zero_one_to_bool, default=False)    # use.mouse.strain
    parser.add_argument('--infer_VJ', type=zero_one_to_bool, default=False)            # infer.VJ
    parser.add_argument('--infer_CDR3', type=zero_one_to_bool, default=False)          # infer.CDR3
    parser.add_argument('--keep_incomplete_chain', type=zero_one_to_bool, default=False)  # keep.incomplete.chain

    # Modeling / run settings
    parser.add_argument('--check_cdr3_mode', type=int, default=1)  # check.cdr3.mode
    parser.add_argument('--start_lg', type=int, default=1)         # start.lg
    parser.add_argument('--end_lg', type=int, default=2)           # end.lg
    parser.add_argument('--renormVJ', type=zero_one_to_bool, default=True)  # logical; R default is NULL
    parser.add_argument('--N_min', type=int, default=10)           # N.min
    parser.add_argument('--build_clones', type=zero_one_to_bool, default=False)  # build.clones
    parser.add_argument('--print_size', type=zero_one_to_bool, default=True)    # print.size

    # Outputs / filenames
    parser.add_argument('--output_stat', type=zero_one_to_bool, default=True)         # output.stat
    parser.add_argument('--output_processed_data', type=zero_one_to_bool, default=False)  # output.processed.data
    parser.add_argument('--filename_output', default="")                              # filename.output -> NULL if ""
    parser.add_argument('--input1_name', default="Input")                             # input1.name
    parser.add_argument('--input2_name', default="")                                  # input2.name -> NULL if ""
    parser.add_argument('--output_format', default="pdf")                             # output.format
    parser.add_argument('--logo_type', default="bits")                                # logo.type

    # Species / model
    parser.add_argument('--species_default', default="HomoSapiens")
    parser.add_argument('--model_default', default="Model_default")

    # Length settings (pass only if provided so R can keep NA)
    parser.add_argument('--set_cdr3a_length', type=int, default=None)
    parser.add_argument('--set_cdr3b_length', type=int, default=None)

    # Verbosity (numeric in R)
    parser.add_argument('--verbose', type=int, default=1)

    args = parser.parse_args()

    # Import R package
    MixTCRviz_pkg = rpackages.importr('MixTCRviz')

    # Build kwargs using exact R parameter names
    r_kwargs = {
        'input1': args.input1,
        'output.path': args.output_path,
        'input2': none_if_empty(args.input2),
        'baseline': none_if_empty(args.baseline_file),

        'chain': args.chain_list_output,
        'interactive.plots': args.interactive_plots,

        'use.allele': args.use_allele,
        'correct.gene.names': args.correct_gene_names,
        'use.mouse.strain': args.use_mouse_strain,

        'check.cdr3.mode': args.check_cdr3_mode,
        'start.lg': args.start_lg,
        'end.lg': args.end_lg,
        'renormVJ': args.renormVJ,
        'N.min': args.N_min,

        'output.stat': args.output_stat,
        'output.processed.data': args.output_processed_data,
        'filename.output': none_if_empty(args.filename_output),

        'plot.title': args.plot_title,
        'set.title': none_if_empty(args.set_title),
        'logo.type': args.logo_type,

        'species.default': args.species_default,
        'model.default': args.model_default,
        'verbose': args.verbose,
        'build.clones': args.build_clones,

        'plot': args.plot,
        'plot.cdr12.motif': args.plot_cdr12_motif,
        'plot.oneline': args.plot_oneline,
        'plot.all.length': args.plot_logo_length,
        'plot.cdr3.norm': args.plot_cdr3_norm,
        'plot.VJ.switch': args.plot_VJ_switch,
        'plot.modelsCombined': args.plot_modelsCombined,

        'label.neg': args.label_neg,
        'label.diag': args.label_diag,
        'plot.sd': args.plot_sd,
        'label.min.fr.input1': args.label_min_fr_input1,
        'label.min.fr.input2': args.label_min_fr_input2,

        'keep.incomplete.chain': args.keep_incomplete_chain,
        'seq.protocol': args.seq_protocol,

        'input1.name': args.input1_name,
        'input2.name': none_if_empty(args.input2_name),
        'output.format': args.output_format,

        'infer.VJ': args.infer_VJ,
        'infer.CDR3': args.infer_CDR3,
        'print.size': args.print_size,

        'ZscoreVJ.thresh': args.ZscoreVJ_thresh,
        'FoldChangeVJ.thresh': args.FoldChangeVJ_thresh,
    }

    # Conditionally include lengths so R can keep default NA if not set
    if args.set_cdr3a_length is not None:
        r_kwargs['set.cdr3a.length'] = args.set_cdr3a_length
    if args.set_cdr3b_length is not None:
        r_kwargs['set.cdr3b.length'] = args.set_cdr3b_length

    # Call R function
# right before:
# MixTCRviz_pkg.MixTCRviz(**r_kwargs)

    r_kwargs = {k: v for k, v in r_kwargs.items() if v is not None}
    MixTCRviz_pkg.MixTCRviz(**r_kwargs)
    ### now convert the rds objects from the stat folder into python pickle dictionaries

    # Define the path to your RDS file
    path_stats = '{0}/stats'.format(args.output_path)
    for rds_file in os.listdir(path_stats):
        if rds_file.split(".")[-1] != 'rds':
            continue
        # Load the RDS file
        data = ro.r['readRDS'](os.path.join(path_stats, rds_file))
        data_python = parse_r_object(data)


        name_pickle_file = os.path.join(path_stats, rds_file.replace('.rds', '.pkl'))
        save_data_to_pickle(name_pickle_file, data_python)               




