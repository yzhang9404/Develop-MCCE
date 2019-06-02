#!/usr/bin/env python
"""
This is an example of logging information. The general rule is to use print() for program output, and logging module for
program progress information.
"""

import logging


def welcome(logger):
    logger.info("Welcome!")
    return


if __name__ == "__main__":
    # Set logging level. DEBUG, INFO, WARNING, ERROR and CRITICAL are in the order from more details to less
    # This has to be declared only once! and before any logging messages!
    logging.basicConfig(level=logging.DEBUG, format='%(name)12s: %(levelname)10s: %(message)s')

    # default is toot logger
    logging.info("Debug message 1")

    # you can use multiple loggers with names
    logger_mc = logging.getLogger("Monte Carlo")
    logger_mc.info("Debug message from Monte Carlo")
    logger_al = logging.getLogger("Analytical")
    logger_al.info("Debug message from Analytical solution")

    # you can send logger to a function
    welcome(logger_mc)