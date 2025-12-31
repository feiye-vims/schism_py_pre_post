import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

def make_session(ignore_proxies=True) -> requests.Session:
    retry = Retry(
        total=6, connect=3, read=3, status=3,
        backoff_factor=1.5,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods={"GET"},
        respect_retry_after_header=True,
    )
    s = requests.Session()
    s.mount("https://", HTTPAdapter(max_retries=retry, pool_connections=32, pool_maxsize=32))
    s.headers.update({"User-Agent": "VIMS-STOFS iv-downloader", "Accept": "text/plain"})
    s.trust_env = not ignore_proxies  # set True only if you *need* system proxy vars
    return s

# Example
with make_session(ignore_proxies=True) as sess:
    r = sess.get(
        "https://waterservices.usgs.gov/nwis/iv",
        params=dict(
            format="rdb", sites="07374000", parameterCd="00060",
            startDT="2018-11-30T00:00:00", endDT="2020-01-03T00:00:00", siteStatus="all",
        ),
        timeout=(5, 180),  # short connect, longer read
    )
    r.raise_for_status()
    print("HTTP", r.status_code, "bytes", len(r.content))
