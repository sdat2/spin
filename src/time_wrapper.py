"""
timeit wrapper.
"""
import time


def timeit(method):
    """
    #  This wrapper is used to make performance comparisons, and is adapted from:
    #  http s://medium.com/@lucasribeiro1990/thank-you-\
    #  so-much-for-this-it-has-been-very-useful-for-me-aeb688ba6f3a
    """

    def timed(*args, **kw):
        ts = time.perf_counter()
        result = method(*args, **kw)
        te = time.perf_counter()
        if "log_time" in kw:
            name = kw.get("log_name", method.__name__.upper())
            kw["log_time"][name] = te - ts
        else:
            print("%r  %2.5f s\n" % (method.__name__, (te - ts)))
        return result

    return timed
