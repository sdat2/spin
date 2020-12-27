import os
import pickle
import src.constants as cst
import src.time_wrapper as twr


@twr.timeit
def print_dict_to_pickle(name="no_name_given", dictionary={}, OUT=cst.DATA_PATH):
    """ Writes """
    if not os.path.exists(OUT):  # make the directory thing
        os.makedirs(OUT)
    with open(OUT + "/" + name + ".pickle", "wb") as handle:
        pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # I can't remember any more why the protocol has to be the highest


@twr.timeit
def read_pickle_to_dict(name="no_name_given", OUT=cst.DATA_PATH):
    """ Reads """
    if not os.path.exists(OUT):  # make the directory thing
        os.makedirs(OUT)
    with open(OUT + "/" + name + ".pickle", "rb") as handle:
        return pickle.load(handle)
    # I can't remember any more why the protocol has to be the highest


@twr.timeit
def circular_write_out(co, sy):
    """ Writes output """
    if not os.path.exists(co.out):  # make the directory thing
        os.makedirs(co.out)
    with open(os.path.join(co.out, co.name + "_co.pickle"), "wb") as handle:
        pickle.dump(co, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # I can't remember any more why the protocol has to be the highest
    with open(os.path.join(co.out, + co.name + "_sy.pickle"), "wb") as handle:
        pickle.dump(sy, handle, protocol=pickle.HIGHEST_PROTOCOL)


@twr.timeit
def read_in(OUTPUT_DIREC, name):
    """load pickles"""
    print("looking for" + os.path.join(OUTPUT_DIREC, name + "_co.pickle"))
    with open(os.path.join(OUTPUT_DIREC, name + "_co.pickle"), "rb") as handle:
        co = pickle.load(handle)
    with open(os.path.join(OUTPUT_DIREC, name + "_sy.pickle"), "rb") as handle:
        sy = pickle.load(handle)
    return co, sy
