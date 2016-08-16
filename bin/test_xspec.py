#! /usr/bin/env python

"""

"""
if __name__ == "__main__":
    import os
    import argparse
    parser = argparse.ArgumentParser(description='see if the import of xspec causes a seg fault')
    parser.add_argument("obsid", type=str, help="Swift XRT obsid to process")
    args=parser.parse_args()
    obsid=args.obsid
    call_xspec()


def call_xspec():
    import xspec