import logging
import shutil
import urllib.request


def get_url(url, fp):
    logging.info("Downloading {0}".format(url))
    with urllib.request.urlopen(url) as resp, open(fp, 'wb') as f:
        shutil.copyfileobj(resp, f)
    return fp


