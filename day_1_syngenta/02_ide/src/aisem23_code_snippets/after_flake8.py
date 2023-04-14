import re
import sys

import pandas as pd


def is_even(number):

    if number % 2 == 0:
        return True
    return False


def get_dates(text):

    return re.findall(r'\d+ \w+ \d{4}', text)


def get_python_versions(text):

    return re.findall(r'Python (\d.+?)<', text)


if __name__ == "__main__":

    print(sys.version)

    my_list = [['a', 'b', 'c'], ['AA', 'BB', 'CC']]

    df = pd.DataFrame(my_list, columns=["data_a", "data_b", "data_c"])

    nums = [1, 2, 3, 4, 5, 6, 7, 8]

    for n in nums:
        if is_even(n):
            print(f"{n} is an even number.")

    # you can scrape the text from the python website using the code below
    # (ensure to import requests module!)
    # taken from: https://www.dataschool.io/web-scraping-with-regex/
    """
    r = requests.get('https://www.python.org/doc/versions/')
    print(r.text[21640:22424])
    """

    text = """
    <h1>Python Documentation by Version</h1>
    <p>Some previous versions of the documentation remain available
    online.  Use the list below to select a version to view.</p>
    <p>For unreleased (in development) documentation, see
    <a class="reference internal" href="#in-development-versions">
    In Development Versions</a>.</p>
    <ul class="simple">
    <li><a class="reference external" href="https://docs.python.org/release/3.11.2/">Python 3.11.2</a>, documentation released on 8 February 2023.</li>
    <li><a class="reference external" href="https://docs.python.org/release/3.11.1/">Python 3.11.1</a>, documentation released on 6 December 2022.</li>
    <li><a class="reference external" href="https://docs.python.org/release/3.11.0/">Python 3.11.0</a>, documentation released on 24 October 2022.</li>
    """

    scraped_dates = get_dates(text)
    print(scraped_dates)

    scraped_versions = get_python_versions(text)
    print(scraped_versions)

    df = pd.DataFrame(zip(scraped_versions, scraped_dates),
                      columns=['Versions', 'Date'])

    print(df.head())
