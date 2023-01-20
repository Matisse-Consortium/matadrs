import argparse
import keyring
import numpy as np
from astroquery.eso import Eso

# Initialize Eso and argparse
parser = argparse.ArgumentParser(description="Query that can be customized to directly load the wanted observations"
											 " to the NoMachine virtual environment.")
eso = Eso()

# Required argument group
requiredName = parser.add_argument_group("required arguments")

# Parser arguments for log in
parser.add_argument("--user", "-u", type=str, required=False, default="MbS", help="changes the default user name")
parser.add_argument("--keyring_removal", "-krm", required=False, default=False,
					help="removes the password stored in keyring.", action="store_true")

# Parser arguments for query
parser.add_argument("--instrument", "-i", type=str, required=False, default="matisse",
					help="changes the default instrument")
parser.add_argument("--target", "-tar", type=str, required=False, default="", help="target source for query")
requiredName.add_argument("--stime", "-st", type=str, required=True,
						  help="start time of observations. Format YYYY-MM-DD")
requiredName.add_argument("--etime", "-et", type=str, required=True, help="end time of observations. Format YYYY-MM-DD")

# Parser argument for help query
parser.add_argument("--help_eso", "-he", required=False, default=False, action="store_true",
					help="shows the query-options for the instrument")

# Parser arguments for action after query
parser.add_argument("--print", "-p", required=False, default=False,
					action="store_true", help="prints the observation table")
parser.add_argument("--print_header", "-ph", required=False, default=False,
					action="store_true", help="prints the headers")
parser.add_argument("--download", "-d", required=False, nargs=2, type=int,
					help="downloads and decompressed the data set."
						 " Takes the data format [<Begin of Download>, <End of Download>]")


# Adds the parser arguments
args = parser.parse_args()

# Log into Eso
eso.login(args.user, store_password=True)
if args.keyring_removal is True:
	keyring.delete_password("astroquery:www.eso.org", args.user)

# Query for the parameters
table = eso.query_instrument(args.instrument, column_filters={"target": args.target, "stime": args.stime,
																"etime": args.etime}, columns=[""], help=args.help_eso)


# Prints the table if print is called
if args.print is True:
	print(table)

# Prints the headers if print_header is called
if args.print_header is True:
	table_headers = eso.get_headers(table["DP.ID"])
	print(table_headers)

data_files = eso.retrieve_data(table['DP.ID'][8:9],
							   destination="/Users/scheuck/Desktop/Programs/PhD/rawData", continuation=True,
							   unzip=False)

