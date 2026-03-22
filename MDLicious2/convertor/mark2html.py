from markdown2 import Markdown
from bs4 import BeautifulSoup

class Mark2HTML:
    def __init__(self, content):
        self.content = content
        self.converter = Markdown(extras=["toc"])

    def convert(self):
        # convert markdown to html
        html = self.converter.convert(self.content)
        
        # modify html with custom behaviour
        html = self.__insert_toc(html)
        html = self.__add_security_to_links(html)
        html = self.__prettify(html)
        
        return html

    def __insert_toc(self, html):
        temp_html = html.split('\n')
        for i in range(len(temp_html)):
            if temp_html[i].find('[toc]') != -1:
                temp_html[i] = html.toc_html
        html = '\n'.join(temp_html)
        return html

    def __add_security_to_links(self, html):
        soup = BeautifulSoup(html, "html.parser")
        for a in soup.find_all("a", href=True):
            href = a["href"]

            # ignore links to current page
            if href.startswith("#"):
                continue

            # ignore links to cfd.university pages
            if not href.startswith(("http://", "https://")):
                continue
            
            # all other pages get get additional security measures
            # to avoid malicious javascript injection (tabnabbing)
            a["target"] = "_blank"
            a["rel"] = "noopener noreferrer"

        return soup.decode(formatter=None)

    def __prettify(self, html):
        soup = BeautifulSoup(html, "html.parser")

        # bs4 removes new lines, adding them back after block level elements
        elements = ["p", "ul", "ol", "h1", "h2", "h3", "h4", "h5", "h6", "iframe", "img", "div", "figure", "table"]
        for tag in soup.find_all(elements):
            tag.insert_after("\n")

        return soup.decode(formatter=None)