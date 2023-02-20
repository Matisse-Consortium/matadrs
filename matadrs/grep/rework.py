import os
import argparse
from astropy.io import fits

# Initialize Eso and argparse
parser = argparse.ArgumentParser(description="Uncompresses the file and shows its data")

# Required argument group
requiredName = parser.add_argument_group("required arguments")

# Required arguments
requiredName.add_argument("--file", "-f", type=str,
                          required=True, help="takes the input in form of a '.fits/.fits.Z'-file")

# Optional arguments
parser.add_argument("--print_info", "-pi", required=False, default=False,
                    action="store_true", help="prints the info of the '.fits'-file")
parser.add_argument("--print_header", "-ph", required=False, default=False,
                    action="store_true", help="prints the header of the '.fits'-file")
# parser.add_argument("--print_data", "-pd", required=False, default=False,
#                    action="store_true", help="prints the data of the '.fits'-file")
parser.add_argument("--print_all_data", "-pad", required=False, default=False,
                    action="store_true", help="prints the data of the '.fits'-file")
parser.add_argument("--print_data", -"pd", type=int, required=False, help="prints the data requested")
args = parser.parse_args()

# Checks if file is in '.Z'-format and if yes uncompresses it
if args.file[-1] == "Z":
    os.system(f'cmd /k "uncompress {args.file}"')
    args.file = args.file[:-2]

# Opens the '.fits' file
hdul = fits.open(args.file, lazy_load_hdus=False)
hdr = hdul[0].header

# Prints info
if args.print_info is True:
    hdul.info()

# Prints headers
if args.print_header is True:
    print(repr(hdr))

# Prints data

# Prints all data
if args.print_all_data is True:
    for i in range(0, 100):
        try:
            print(hdul[i].data)
            print()
        except:
            pass
