#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/sequana/sequana
#  documentation: http://sequana.readthedocs.io
#
##############################################################################
""".. rubric:: misc utilities"""
import asyncio

import aiohttp
from tqdm.asyncio import tqdm

import colorlog

from sequana.lazy import numpy as np

logger = colorlog.getLogger(__name__)


__all__ = ["textwrap", "wget", "download", "findpos", "normpdf", "multiple_downloads"]


def normpdf(x, mu, sigma):
    """Return the normal pdf evaluated at *x*; args provides *mu*, *sigma*"

    .. note:: same as scipy.stats.norm but implemented to avoid scipy dependency
    """

    return 1.0 / (np.sqrt(2 * np.pi) * sigma) * np.exp(-0.5 * (1.0 / sigma * (x - mu)) ** 2)


def textwrap(text, width=80, indent=0):
    """Wrap a string with 80 characters

    :param text: input text
    :param width: (defaults to 80 characters)
    :param indent: possible indentation (0 by default)

    """
    if indent == 0:
        indent = ""
    else:
        indent = " " * indent
    data = [indent + text[i * width : (i + 1) * width :] for i in range(len(text) // width + 1)]
    return "\n".join(data)


def wget(link, output):
    """Retrieve a file from internet.

    :param str link: a valid URL
    :param str output: the output filename

    .. warning:: no sanity check of any kind for now
    """
    try:
        from urllib import urlretrieve
    except:
        from urllib.request import urlretrieve
    urlretrieve(link, filename=output)


def findpos(seq, chr):
    """Find position(s) of a substring into a longer string.

    Note that this function is a generator::

        >>> list(findpos("AACCGGAAGGTT", "GG"))
        [4,8]

    """
    N = len(chr)
    for i, dummy in enumerate(seq):
        if seq[i : i + N] == chr:
            yield i


def download(url, output):
    """Download a file from a given URL using asynchronous HTTP requests.

    :param str url: The URL from which to download the file.
    :param str output: The path and filename where the downloaded file will be saved.

    Raises a KeyboardInterrupt or asyncio.TimeoutError if the download is interrupted or takes too long.
    In such cases, it logs a message, removes partially downloaded files, and logs a critical message.
    """


    files_to_download = [(url, output, 0)]
    try:  # try an asynchrone downloads
        multiple_downloads(files_to_download)
    except (KeyboardInterrupt, asyncio.TimeoutError): # pragma: no cover
        logger.info("The download was interrupted or network was too slow. Removing partially downloaded files")
        for values in files_to_download:
            filename = values[1]
            Path(filename).unlink()
        logger.critical(
            "Keep going but your pipeline will probably not be fully executable since images could not be downloaded"
        )



# copied from sequana_pipetools
def multiple_downloads(files_to_download, timeout=3600):
    async def download(session, url, name, position):
        async with session.get(url, timeout=timeout) as resp:
            with tqdm.wrapattr(
                open(name, "wb"),
                "write",
                miniters=1,
                desc=url.split("/")[-1],
                total=int(resp.headers.get("content-length", 0)),
                position=position,
            ) as fout:
                async for chunk in resp.content.iter_chunked(4096):
                    fout.write(chunk)

    async def download_all(files_to_download):
        """data_to_download is a list of tuples
        each tuple contain the url to download, its output name, and a unique
        position for the progress bar."""
        async with aiohttp.ClientSession(connector=aiohttp.TCPConnector(limit=10)) as session:
            await asyncio.gather(*(download(session, *data) for data in files_to_download))

    asyncio.run(download_all(files_to_download))
