# Place all your constants here
import os

# Note: constants should be UPPER_CASE
constants_path = os.path.realpath(__file__)
SRC_PATH = os.path.dirname(constants_path)
PROJECT_PATH = os.path.dirname(SRC_PATH)
DATA_PATH = os.path.join(PROJECT_PATH, "io")
NC_PATH = os.path.join(DATA_PATH, "nc")
