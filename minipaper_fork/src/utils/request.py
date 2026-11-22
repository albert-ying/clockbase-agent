import os
import hashlib
import asyncio
import aiohttp
import shutil
import ssl
import base64
import requests
from aiohttp_socks import ProxyConnector, ProxyType
from functools import wraps
from typing import Union
from playwright.async_api import async_playwright, ProxySettings
from playwright.sync_api import sync_playwright
import time
from datetime import datetime
from typing import Optional, List
from pathlib import Path
from dotenv import load_dotenv

# local module
from utils.logger import logger
from utils.helpers import retry as sync_retry
from utils.http_client import safe_get, safe_post, get_proxies  # timeout + retry + proxy fix


load_dotenv()


def async_retry(max_retries:int=5):
    """Retry decorator for async functions.
    Args:
        max_retries (int): The maximum number of retries.
    """
    def retry_decorator(func):
        @wraps(func)
        async def wrapper(*args, **kwargs):
            for i in range(max_retries):
                try:
                    return await func(*args, **kwargs)
                except Exception as e:
                    logger.error(f'Error in {func.__name__}: {e}')
                    # raise e
                    await asyncio.sleep(30 * (i + 1))
            raise Exception(f'Failed to execute {func.__name__} after {max_retries} retries.')
        return wrapper
    return retry_decorator


async def fetch(
        url:str, 
        data:Union[dict, aiohttp.FormData] = None, 
        cls_name:str = None, 
        timeout:float=1800, 
        session_kw:dict=None, 
        request_kw:dict=None,
        semaphore:asyncio.Semaphore=None,
        return_type:str='json',
        method:str='GET'
    ):
        """This function is used to fetch data from a given url. Support GET and POST methods.
        If the method is GET, the data will be sent as a query string.
        If the method is POST, the data will be sent as a JSON object.
        If the method is POST and the data is not a JSON object, the data will be sent as a form data.
        Args:
            url (str): The url to fetch data from.
            data (dict): The data to send to the url.
            cls_name (str): The name of the class that is making the request.
            timeout (float): The timeout for the request.
            session_kw (dict): The keyword arguments for the session.
            request_kw (dict): The keyword arguments for the request.
            semaphore (asyncio.Semaphore): The semaphore for the request.
            return_type (str): The type of the return value.
            method (str): The method of the request.
        """
        data = data or dict()
        cls_name = cls_name or 'fetch'
        session_kw = session_kw or dict()
        request_kw = request_kw or dict()
        ssl_context = ssl.create_default_context()
        ssl_context.check_hostname = False
        ssl_context.verify_mode = ssl.CERT_NONE
        if 'connector' in session_kw:
            connector = session_kw.pop('connector')
        else:
            connector = aiohttp.TCPConnector(ssl=ssl_context)
        timeout = aiohttp.ClientTimeout(total=timeout)
        semaphore = semaphore or asyncio.Semaphore(1)
        async with semaphore:
            async with aiohttp.ClientSession(connector=connector, timeout=timeout, **session_kw) as session:
                try:                    
                    if method == 'GET':
                        request_kw.update({
                            'url': url,
                            'params': data
                        })
                        async with session.get(**request_kw) as resp:
                            if resp.status == 200:
                                if return_type == 'json':
                                    result = await resp.json()
                                elif return_type == 'text':
                                    result = await resp.text()
                                elif return_type == 'content':
                                    result = await resp.read()
                            else:
                                raise ValueError(f'{cls_name} Request Failed: {resp.status}.')
                            return result
                    else:
                        if isinstance(data, dict):
                            request_kw.update({
                                'url': url,
                                'json': data
                            })
                        else:
                            request_kw.update({
                                'url': url,
                                'data': data
                            })
                        async with session.post(**request_kw) as resp:
                            if resp.status == 200:
                                if return_type == 'json':
                                    result = await resp.json()
                                elif return_type == 'text':
                                    result = await resp.text()
                                elif return_type == 'content':
                                    result = await resp.read()
                            else:
                                raise ValueError(f'{cls_name} Request Failed: {resp.status}.')
                            return result
                except asyncio.TimeoutError as e:
                    logger.info(f'{cls_name} Request Timeout.')
                    raise e
                except aiohttp.ClientError as e:
                    logger.info(f"{cls_name} Request Error.")
                    raise e


class AsyncPDFDownloader:
    def __init__(self, download_dir: str = "downloads", headless: bool = True):
        """
        Initializes the PDF downloader.
        :param download_dir: The directory to save PDFs.
        """
        self.download_dir = Path(download_dir)
        if not self.download_dir.exists():
            self.download_dir.mkdir(parents=True, exist_ok=True)
        print(f"Download directory set to: {self.download_dir.resolve()}")
        self.playwright = None
        self.browser = None
        self.headless = headless

    async def __aenter__(self):
        self.playwright = await async_playwright().start()
        self.browser = await self.playwright.chromium.launch(
            headless=self.headless,
            args=[
                "--disable-web-security",
                "--disable-features=IsolateOrigins,site-per-process",
            ]
        )
        return self

    async def __aexit__(self, exc_type, exc_val, exc_tb):
        if self.browser:
            await self.browser.close()
        if self.playwright:
            await self.playwright.stop()

    def _generate_filename(self, url: str) -> str:
        """
        Generates a unique filename based on the URL and timestamp.
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        url_hash = hashlib.md5(url.encode()).hexdigest()[:8]
        return f"pdf_{timestamp}_{url_hash}.pdf"

    async def download_pdf(self, url: str, filename: Optional[str] = None) -> Optional[str]:
        """
        Asynchronously downloads a single PDF file.
        """
        if not self.browser:
            raise RuntimeError("Browser not initialized. Use this class with an 'async with' statement.")

        if not filename:
            filename = self._generate_filename(url)

        filepath = self.download_dir / filename
        print(f"\nStarting PDF download...")
        print(f"Target URL: {url}")
        print(f"Expected save location: {filepath}")

        context = None
        page = None
        try:
            print("\nCreating new browser context...")
            context = await self.browser.new_context(
                accept_downloads=True,
                user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
            )

            print("\nCreating new page...")
            page = await context.new_page()

            print("\nSetting browser headers...")
            await page.set_extra_http_headers(
                {
                    "Accept": "application/pdf,application/x-pdf,*/*",
                    "Accept-Language": "en-US,en;q=0.9",
                    "Sec-Fetch-Dest": "document",
                    "Sec-Fetch-Mode": "navigate",
                    "Sec-Fetch-Site": "none",
                    "Sec-Fetch-User": "?1",
                    "Upgrade-Insecure-Requests": "1",
                }
            )

            print("\nPreparing to download...")
            max_retries = 3
            retry_count = 0

            while retry_count < max_retries:
                try:
                    async with page.expect_download(timeout=30000) as download_info:
                        try:
                            print(f"Navigating to URL: {url}")
                            response = await page.goto(
                                url,
                                wait_until="commit",
                                timeout=30000
                            )
                            await page.wait_for_load_state("networkidle")
                            
                            if response:
                                print(f"Response status: {response.status}")
                                print(f"Response content type: {response.headers.get('content-type', 'unknown')}")
                                ...
                        except Exception as e:
                            print(f"Navigation exception: {e}")
                            # Continue execution as the download may have already started

                    print("Waiting for download to start...")
                    download = await download_info.value
                    print(f"Suggested filename: {download.suggested_filename}")

                    print("Waiting for download to complete...")
                    path = await download.path()
                    print(f"Download temporary location: {path}")

                    filepath.parent.mkdir(parents=True, exist_ok=True)

                    print(f"Moving file to target location: {filepath}")
                    shutil.move(path, filepath)

                    print("\n=== PDF downloaded successfully ===")
                    return str(filepath)

                except Exception as e:
                    retry_count += 1
                    print(f"\nDownload attempt {retry_count} failed: {str(e)}")
                    if retry_count < max_retries:
                        print(f"Retrying after 3 seconds...")
                        await asyncio.sleep(3)
                    else:
                        print("Maximum retries reached.")
                        return None

        except Exception as e:
            print(f"\n!!! An unexpected error occurred !!!")
            print(f"Error type: {type(e).__name__}")
            print(f"Error message: {str(e)}")
            return None

        finally:
            if page:
                await page.close()
            if context:
                await context.close()

    async def download_multiple_pdfs(self, urls: List[str], filenames: Optional[List[str]] = None) -> List[Optional[str]]:
        """
        Asynchronously downloads multiple PDF files.
        :param urls: List of PDF URLs.
        :param filenames: Optional list of filenames.
        :return: List of downloaded file paths.
        """
        if filenames is None:
            filenames = [None] * len(urls)
        elif len(filenames) != len(urls):
            raise ValueError("The number of URLs and filenames must be the same.")

        # Create a list of download tasks
        tasks = [
            self.download_pdf(url, filename)
            for url, filename in zip(urls, filenames)
        ]

        # Run all download tasks concurrently
        results = await asyncio.gather(*tasks)
        return results

    def verify_pdf(self, filepath: str) -> bool:
        """
        Verifies if the downloaded file is a valid PDF.
        :param filepath: Path to the PDF file.
        :return: True if the file is a valid PDF, False otherwise.
        """
        print(f"\nVerifying PDF file: {filepath}")

        path = Path(filepath)
        if not path.exists():
            print("File does not exist!")
            return False

        # Check file size
        file_size = path.stat().st_size
        print(f"File size: {file_size} bytes")

        if file_size < 100:
            print("File is too small, likely invalid.")
            return False

        # Check if the file starts with the PDF magic number
        try:
            with open(filepath, "rb") as f:
                header = f.read(4)
                is_pdf = header == b"%PDF"
                print(f"PDF header check: {'Passed' if is_pdf else 'Failed'}")
                return is_pdf
        except Exception as e:
            print(f"An error occurred while checking the file: {str(e)}")
            return False


class PDFDownloader:
    def __init__(self, download_dir: str = "downloads", headless: bool = True):
        """
        Initializes the PDF downloader.
        :param download_dir: The directory to save PDFs.
        """
        self.download_dir = Path(download_dir)
        if not self.download_dir.exists():
            self.download_dir.mkdir(parents=True, exist_ok=True)
        
        self.playwright = sync_playwright().start()
        self.browser = self.playwright.chromium.launch(
            headless=headless,
            args=[
                "--disable-web-security",
                "--disable-features=IsolateOrigins,site-per-process",
            ]
        )
        print(f"Download directory set to: {self.download_dir.resolve()}")

    def close(self):
        if self.browser:
            self.browser.close()
        if self.playwright:
            self.playwright.stop()

    def _generate_filename(self, url: str) -> str:
        """
        Generates a unique filename based on the URL and timestamp.
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        url_hash = hashlib.md5(url.encode()).hexdigest()[:8]
        return f"pdf_{timestamp}_{url_hash}.pdf"

    def download_pdf(self, url: str, filename: Optional[str] = None) -> Optional[str]:
        """
        Synchronously downloads a single PDF file.
        """
        if not filename:
            filename = self._generate_filename(url)

        filepath = self.download_dir / filename
        print(f"\nStarting PDF download...")
        print(f"Target URL: {url}")
        print(f"Expected save location: {filepath}")

        context = None
        page = None
        try:
            print("\nCreating new browser context...")
            context = self.browser.new_context(
                accept_downloads=True,
                user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
            )

            print("\nCreating new page...")
            page = context.new_page()

            print("\nSetting browser headers...")
            page.set_extra_http_headers(
                {
                    "Accept": "application/pdf,application/x-pdf,*/*",
                    "Accept-Language": "en-US,en;q=0.9",
                    "Sec-Fetch-Dest": "document",
                    "Sec-Fetch-Mode": "navigate",
                    "Sec-Fetch-Site": "none",
                    "Sec-Fetch-User": "?1",
                    "Upgrade-Insecure-Requests": "1",
                }
            )

            print("\nPreparing to download...")
            max_retries = 3
            retry_count = 0

            while retry_count < max_retries:
                try:
                    with page.expect_download(timeout=30000) as download_info:
                        try:
                            print(f"Navigating to URL: {url}")
                            response = page.goto(
                                url,
                                wait_until="commit",
                                timeout=30000
                            )
                            
                            if response:
                                print(f"Response status: {response.status}")
                                print(f"Response content type: {response.headers.get('content-type', 'unknown')}")
                        except Exception as e:
                            print(f"Navigation exception: {e}")

                    print("Waiting for download to start...")
                    download = download_info.value
                    print(f"Suggested filename: {download.suggested_filename}")

                    print("Waiting for download to complete...")
                    path = download.path()
                    print(f"Download temporary location: {path}")

                    filepath.parent.mkdir(parents=True, exist_ok=True)

                    print(f"Moving file to target location: {filepath}")
                    shutil.move(path, filepath)

                    print("\n=== PDF downloaded successfully ===")
                    return str(filepath)

                except Exception as e:
                    retry_count += 1
                    print(f"\nDownload attempt {retry_count} failed: {str(e)}")
                    if retry_count < max_retries:
                        print(f"Retrying after 3 seconds...")
                        time.sleep(3)
                    else:
                        print("Maximum retries reached.")
                        return None

        except Exception as e:
            print(f"\n!!! An unexpected error occurred !!!")
            print(f"Error type: {type(e).__name__}")
            print(f"Error message: {str(e)}")
            return None

        finally:
            if page:
                page.close()
            if context:
                context.close()

def download_pdf_after_second_visit(url, output_filename="downloaded_pdf.pdf"):
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=False, args=['--disable-web-security', '--disable-features=IsolateOrigins,site-per-process']) # Set headless=False to observe browser behavior
        context = browser.new_context(
            accept_downloads=True,
            user_agent="Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
        )   
        page = context.new_page()

        try:
            print(f"Performing first visit to complete challenge and set cookies: {url}")
            # First visit: navigate to URL and wait for network idle to ensure JavaScript (including POW) execution completes
            page.goto(url, wait_until='networkidle', timeout=90000) # Increased timeout for complex JS loading
            print("First visit complete. Waiting to ensure all background tasks finish...")
            page.wait_for_timeout(5000) # Extra 5-second wait to ensure cookies and other state are stable

            # Check for any redirects or errors after the first visit
            if "about:blank" in page.url or "Preparing to download" in page.title():
                print("Page is still preparing download or showing blank after first visit. Proceeding with second visit.")
            else:
                print(f"Page URL after first visit: {page.url}")
                print(f"Page title after first visit: {page.title()}")
                # If PDF is already displayed on the first visit, handle it directly
                try:
                    embed_element = page.wait_for_selector('embed[type="application/pdf"]', state='attached', timeout=5000)
                    current_src = embed_element.get_attribute('src')
                    if current_src and current_src != 'about:blank':
                        print("PDF already loaded on first visit, handling directly.")
                        return _handle_pdf_src(page, current_src, output_filename)
                except Exception:
                    # PDF not found or src not updated, proceed with second visit
                    pass

            print(f"Performing second visit to retrieve PDF: {url}")
            # Second visit: navigate to the same URL again, browser will carry cookies set during the first visit
            with page.expect_download() as download_info:
                page.goto(url, wait_until='domcontentloaded', timeout=90000) # Second load, wait for DOM content
                pdf_src = None
                start_time = time.time()
                timeout = 5000 # Longer timeout to wait for src attribute update

                print(f"Waiting for embed element src attribute to update (max {timeout} seconds) after second visit...")
                # while (time.time() - start_time) < timeout:
                #     try:
                #         locator = page.locator('embed[type="application/x-google-chrome-pdf"]')
                #         if locator:
                #             current_src = locator.get_attribute('src')
                #             if current_src and current_src != 'about:blank':
                #                 pdf_src = current_src
                #                 print(f"Detected new src: {pdf_src}")
                #                 break
                #         else:
                #             print("Embed element not found, continuing to wait...")
                #     except Exception as e:
                #         print(f"Error while getting src: {e}")
                #     time.sleep(1) # Check every second
                
                # if not pdf_src:
                #     print(f"Failed to detect embed element src attribute update within {timeout} seconds.")
                #     print(f"Current page URL (second visit): {page.url}")
                #     print(f"Current page title (second visit): {page.title()}")
                #     return False

                body_html = page.locator('body').inner_html()
                print(body_html)

                download = download_info.value
                if download.suggested_filename == output_filename:
                    print(f"PDF successfully downloaded to {output_filename}")
                    return True
                else:
                    print(f"Downloaded filename is incorrect: {download.suggested_filename}")


            return _handle_pdf_src(page, pdf_src, output_filename)

        except Exception as e:
            print(f"Error occurred while downloading PDF: {e}")
            return False
        finally:
            browser.close()

def _handle_pdf_src(page, pdf_src, output_filename):
    """Helper function to handle PDF download based on src type."""
    if pdf_src.startswith('data:application/pdf;base64,'):
        # Handle Data URI
        print("Detected PDF as Data URI.")
        try:
            encoded_data = pdf_src.split(',', 1)[1] # split(..., 1) ensures splitting only at the first comma
            pdf_content = base64.b64decode(encoded_data)
            with open(output_filename, "wb") as f:
                f.write(pdf_content)
            print(f"PDF (Data URI) successfully downloaded to {output_filename}, size: {len(pdf_content)} bytes.")
            return True
        except Exception as e:
            print(f"Error decoding Data URI or saving file: {e}")
            return False
    elif pdf_src.startswith('http://') or pdf_src.startswith('https://'):
        # Handle direct URL
        print(f"Detected PDF as a direct URL: {pdf_src}")
        try:
            response = page.request.get(pdf_src)
            if response.ok:
                pdf_content = response.body()
                with open(output_filename, "wb") as f:
                    f.write(pdf_content)
                print(f"PDF successfully downloaded to {output_filename}, size: {len(pdf_content)} bytes.")
                return True
            else:
                print(f"Failed to download PDF from {pdf_src}, status code: {response.status}")
                return False
        except Exception as e:
            print(f"Error downloading with page.request.get: {e}")
            return False
    else:
        print(f"Unrecognized PDF source type: {pdf_src}")
        return False


class SemanticScholarRequester:
    def __init__(self, task_queue:asyncio.Queue, semaphore:asyncio.Semaphore=None):
        self.session = None
        self.semaphore = semaphore or asyncio.Semaphore(1)
        self.api_url = 'https://api.semanticscholar.org/graph/v1/paper/search'
        self.task_queue = task_queue

    @async_retry(max_retries=10)
    async def get_papers(self, query:str, start_date:datetime, end_date:datetime, offset:int=0, limit:int=100):
        params = {
            'query': query,
            'openAccessPdf': 'true',
            'publicationDateOrYear': f'{start_date.strftime("%Y-%m-%d")}:{end_date.strftime("%Y-%m-%d")}',
            'fields': 'abstract,authors,citationCount,citationStyles,corpusId,externalIds,fieldsOfStudy,influentialCitationCount,isOpenAccess,journal,openAccessPdf,paperId,publicationDate,publicationTypes,publicationVenue,referenceCount,s2FieldsOfStudy,title,url,venue,year',
            'offset': offset,
            'limit': limit
        }
        request_kw = {
            'headers': {
                'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36',
                'Accept': 'application/json',
                'Accept-Language': 'en-US,en;q=0.9',
                'Accept-Encoding': 'gzip, deflate, br',
                'Connection': 'keep-alive',
                'Referer': 'https://www.semanticscholar.org/',
                'Origin': 'https://www.semanticscholar.org',
            },
            'proxy': f'http://{os.getenv("PROXY_USER")}:{os.getenv("PROXY_PWD")}@{os.getenv("PROXY_HOST")}:{os.getenv("PROXY_PORT")}'
        }
        response = await fetch(self.api_url, params, 'SemanticScholarRequester', semaphore=self.semaphore, method='GET', request_kw=request_kw)
        logger.info(f'Requesting papers from {start_date.strftime("%Y-%m-%d")} to {end_date.strftime("%Y-%m-%d")}')
        return response


    async def loop_get_papers(self):
        while True:
            await asyncio.sleep(1)
            paper = await self.task_queue.get()
            await self.get_papers(paper['query'], paper['start_date'], paper['end_date'])
            if self.task_queue.empty():
                break


    @sync_retry(times=10, delay=30)
    def get_pdf_url_by_doi(self, doi:str):
        api_url = f'https://api.semanticscholar.org/graph/v1/paper/{doi}'


        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0.0.0 Safari/537.36',
            'Accept': 'application/json',
            'Accept-Language': 'en-US,en;q=0.9',
            'Accept-Encoding': 'gzip, deflate, br',
            'Connection': 'keep-alive',
            'Referer': 'https://www.semanticscholar.org/',
            'Origin': 'https://www.semanticscholar.org'
        }
        proxies = get_proxies()  # correctly-keyed http/https dict or None
        params = {'fields': 'openAccessPdf'}
        response = safe_get(api_url, headers=headers, proxies=proxies, params=params)
        if response.status_code != 200:
            return None
        pdf_url = response.json()['openAccessPdf'].get('url')
        if pdf_url:
            return pdf_url
        else:
            return None



if __name__ == "__main__":
    # pmc = 'PMC9705299'
    # link_tag = 'pdf/41467_2022_Article_34854.pdf'
    # url = f'https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/{link_tag}'
    # headers = {
    #     'authority': 'pmc.ncbi.nlm.nih.gov',
    #     'method': 'GET',
    #     'path': f'/articles/{pmc}/{link_tag}',
    #     "referer": f"https://pmc.ncbi.nlm.nih.gov/articles/{pmc}",
    #     'scheme': 'https',
    #     'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/137.0.0.0 Safari/537.36',
    #     "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7",
    #     "Accept-Encoding": "gzip, deflate, br, zstd",
    #     "Accept-Language": "en-US,en;q=0.9",
    #     "sec-ch-ua-platform": "Windows",
    #     "Sec-Fetch-Dest": "document",
    #     "Sec-Fetch-Mode": "navigate",
    #     "Sec-Fetch-Site": "same-origin",
    #     "Sec-Fetch-User": "?1",
    #     "Upgrade-Insecure-Requests": "1"
    # }


    # url = 'https://pmc.ncbi.nlm.nih.gov/articles/PMC9705299/pdf/41467_2022_Article_34854.pdf'
    # download_pmc_pdf(url, output_filename="downloaded_pdf.pdf", headless=False, headers=headers)
    # download_pdf_after_second_visit(url, "PMC9705299_double_visit.pdf")

    # requester = SemanticScholarRequester(asyncio.Queue())
    # pdf_url =  requester.get_pdf_url_by_doi('10.1038/s41467-022-34854-w')
    # print(pdf_url)

    downloader = PDFDownloader()
    downloader.download_pdf('https://pmc.ncbi.nlm.nih.gov/articles/PMC9705299/pdf/41467_2022_Article_34854.pdf')




