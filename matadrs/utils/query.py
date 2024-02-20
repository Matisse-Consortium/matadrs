import getpass
from pathlib import Path
from typing import Optional, List

import cgi
import keyring
import requests
from datetime import datetime, timedelta
from astropy.table import vstack
from astroquery.eso import Eso


def query_archive(
        user_name: str, instrument: str,
        target: str, nights: List[str],
        columns: Optional[List[str]] = [""],
        remove_password: Optional[bool] = False,
        row_limit: Optional[int] = 1500,
        help: Optional[bool] = False,
        download: Optional[bool] = False,
        save_dir: Optional[str] = None
        ) -> None:
    """Queries the ESO's raw data archive.

    Parameters
    ----------
    user_name : str
        The user name for the ESO archive.
    instrument : str
        The instrument to query the archive for.
    target : str
        The target to query the archive for.
    nights : list of str
        The nights to be queried.
    row_limit : int, optional
        The maximum number of rows to return.
        If -1, all rows are returned (takes quite some time).
    remove_password : bool, optional
        If True, the password will be removed from the keyring.
    help : bool, optional
        If True, the help message will be printed.
    download : bool, optional
        If True, the data will be downloaded.
    save_dir : str, optional
        The directory to save the data to. The data will be sorted in
        the folders pertaining to their array configuration and night.
    """
    eso, server = Eso(), "astroquery:www.eso.org"
    eso.ROW_LIMIT = row_limit
    eso.login(user_name, store_password=True)
    password = keyring.get_password(server, user_name)
    if password is None:
        password = getpass.getpass(f"Reenter password: ")

    if remove_password:
        keyring.delete_password(server, user_name)

    for night in nights:
        next_night = (datetime.strptime(night, "%Y-%m-%d")+timedelta(days=1)).strftime("%Y-%m-%d")
        column_filters = {"target": target, "stime": night, "etime": next_night}
        sci = eso.query_instrument(instrument, column_filters=column_filters,
                                   columns=columns, help=help)
        sci_headers = eso.get_headers(sci['DP.ID'])
        sci_starts = np.unique(sci_headers["HIERARCH ESO TPL START"]).data
        sci_name = sci_headers["HIERARCH ESO OBS TARG NAME"][0]
        array_configuration = sci_headers["HIERARCH ESO ISS BASELINE"][0]

        del column_filters["target"]
        night = eso.query_instrument(instrument, column_filters=column_filters,
                                     columns=columns, help=help)
        night_headers = eso.get_headers(night['DP.ID'])
        tpl_starts = np.unique(night_headers["HIERARCH ESO TPL START"]).data
        ind = int(np.where(tpl_starts == sci_starts)[0].squeeze())

        calibrators = []
        for i in [ind-1, ind+1]:
            obs = night_headers[night_headers["HIERARCH ESO TPL START"] == tpl_starts[i]]
            dpr_type = obs["HIERARCH ESO DPR CATG"][0]
            if dpr_type == "SCIENCE":
                if obs["HIERARCH ESO OBS TARG NAME"][0] == sci_name:
                    print("[WARNING]: Sequential observation of the same SCI.")
                continue
            else:
                ob_name = obs["HIERARCH ESO OBS NAME"][0]
                calibrators.append(obs)
                if sci_name.replace("_", "").lower() in ob_name.replace("_", "").lower():
                    if "LN" in ob_name.upper():
                        break
    
        if download:
            download_files = vstack([sci_headers, *calibrators])["DP.ID"]
            destination = save_dir / array_configuration / night

            retrieve_data(download_files, destination=destination, with_calib="raw",
                          username=user_name, password=password)


def getToken(username, password):
    """Token based authentication to ESO: provide username and password to receive back a JSON Web Token."""
    if username==None or password==None:
        return None
    token_url = "https://www.eso.org/sso/oidc/token"
    token = None
    try:
        response = requests.get(token_url,
                            params={"response_type": "id_token token", "grant_type": "password",
                                    "client_id": "clientid",
                                    "username": username, "password": password})
        token_response = json.loads(response.content)
        token = token_response['id_token']+'=='
    except NameError as e:
        print(e)
    except:
        print("*** AUTHENTICATION ERROR: Invalid credentials provided for username %s" %(username))
    
    return token


def getDispositionFilename( response ):
    """Get the filename from the Content-Disposition in the response's http header"""
    contentdisposition = response.headers.get('Content-Disposition')
    if contentdisposition == None:
        return None
    _, params = cgi.parse_header(contentdisposition)
    filename = params["filename"]
    return filename


def writeFile(response, destination):
    """Write on disk the retrieved file"""
    if response.status_code == 200:
        filename = getDispositionFilename(response)
        with open(destination / filename, 'wb') as f:
            f.write(response.content)
        return filename 


def retrieve_data(download_files, destination, with_calib,
                  username, password):
    for name in download_files["DP.ID"]:
        file_url = f'https://dataportal.eso.org/dataportal_new/file/{name}'


        username = input("Type your ESO username: ")
        password=getpass.getpass(prompt="%s user's password: "%(username), stream=None)
        token = getToken(username, password)

        headers = None
        if token!=None:
            headers = {"Authorization": "Bearer " + token}
            response = requests.get(file_url, headers=headers)
            filename = writeFile(response, destination)
            if filename:
                print("Saved file: %s" % (destination / filename))
            else:
                print("Could not get file (status: %d)" % (response.status_code))
        else:
            print("Could not authenticate")



if __name__ == "__main__":
    query_archive("MbS", instrument="matisse",
                  target="hd142527", nights=[""],
                  save_dir=Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/matisse/GTO/hd142527/raw"))
