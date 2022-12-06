import os
import sys
import shutil
import importlib.util
import subprocess

from setuptools import setup, find_packages


# NOTE: Avoids errors during pre-install collection of wxPython
# (which is required for mat-tools)
if importlib.util.find_spec("attrdict") is None:
    subprocess.run([sys.executable, "-m", "pip", "install", "attrdict==2.0.1"],
                   check=True)

setup(
    name='matadrs',
    version='0.1',
    packages=find_packages(include=["matadrs"]),
    package_dir={"": "src"},
    python_requires=">=3.8, <=3.10",
    install_requires=[
        # Requirements for mat_tools
        "astropy==5.1.1",
        "astroquery==0.4.6",
        "attrdict==2.0.1",
        "beautifulsoup4==4.11.1",
        "certifi==2022.9.24",
        "charset-normalizer==2.1.1",
        "contourpy==1.0.6",
        "cycler==0.11.0",
        "docutils==0.19",
        "et-xmlfile==1.1.0",
        "fonttools==4.38.0",
        "html5lib==1.1",
        "idna==3.4",
        "importlib-metadata==5.1.0",
        "jaraco.classes==3.2.3",
        "keyring==23.11.0",
        "kiwisolver==1.4.4",
        "matadrs==0.1",
        "matplotlib==3.6.2",
        "more-itertools==9.0.0",
        "numpy==1.23.5",
        "ObjectListView==1.3.1",
        "openpyxl==3.0.10",
        "packaging==21.3",
        "Pillow==9.3.0",
        "pyerfa==2.0.0.1",
        "pyparsing==3.0.9",
        "python-dateutil==2.8.2",
        "pyvo==1.4",
        "PyYAML==6.0",
        "requests==2.28.1",
        "Shapely==1.8.5.post1",
        "six==1.16.0",
        "skycalc-cli==1.4",
        "soupsieve==2.3.2.post1",
        "statistics==1.0.3.5",
        "tqdm==4.64.1",
        "urllib3==1.26.13",
        "webencodings==0.5.1",
        "wxPython==4.2.0",
        "zipp==3.11.0",

        # Direct requirements for matadrs
        "wget==3.2",

    ]
)

if importlib.util.find_spec("mat_tools") is None:
    import wget

    mat_tools_remote = "https://gitlab.oca.eu/MATISSE/tools/-/archive/master/tools-master.zip"
    reduction_dir = os.path.abspath("src/matadrs/libs/reduction")
    mat_tools_dir = os.path.join(reduction_dir, "tools-master")
    zip_file_path = mat_tools_dir+".zip"

    if not os.path.exists(zip_file_path):
        wget.download(mat_tools_remote, zip_file_path)
    if not os.path.exists(mat_tools_dir):
        subprocess.run(["tar", "-xf", zip_file_path], cwd=reduction_dir, check=True)

    subprocess.run([sys.executable, "-m", "pip", "-e", "install", "."], cwd=mat_tools_dir,
                   check=True)
    shutil.rm(zip_file_path)

