import os
import pickle


def print_dict_to_pickle(name="no_name_given", dictionary={}, OUT="Pickle_Dump"):
    """ Writes """
    if not os.path.exists(OUT):  # make the directory thing
        os.makedirs(OUT)
    with open(OUT + "/" + name + ".pickle", "wb") as handle:
        pickle.dump(dictionary, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # I can't remember any more why the protocol has to be the highest


def read_pickle_to_dict(name="no_name_given", OUT="Pickle_Dump"):
    """ Reads """
    if not os.path.exists(OUT):  # make the directory thing
        os.makedirs(OUT)
    with open(OUT + "/" + name + ".pickle", "rb") as handle:
        return pickle.load(handle)
    # I can't remember any more why the protocol has to be the highest


def circular_write_out(co, sy):
    """Writes output"""
    if not os.path.exists(co.out):  # make the directory thing
        os.makedirs(co.out)
    with open(co.out + "/" + co.name + "_co.pickle", "wb") as handle:
        pickle.dump(co, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # I can't remember any more why the protocol has to be the highest
    with open(co.out + "/" + co.name + "_sy.pickle", "wb") as handle:
        pickle.dump(sy, handle, protocol=pickle.HIGHEST_PROTOCOL)
