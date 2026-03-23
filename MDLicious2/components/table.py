from markdown2 import Markdown
from bs4 import BeautifulSoup

from MDLicious2.components.base import Component
from MDLicious2.convertor.captionExtractor import ComponentType

class Table(Component):
    def __init__(self, content):
        super().__init__(content)

    def match(self, index):
        line = self.content[index]
        is_table = line.find('<!-- table') != -1
        return is_table

    def convert(self, index):
        # check how many lines the table has
        is_table = True
        num_lines = 0
        while is_table:
            if index + num_lines >= len(self.content):
                is_table = False
                break
            
            line = self.content[index + num_lines]
            
            if line.find('<!-- table') != -1:
                num_lines += 1
                continue
            
            if len(line.strip()) > 0:
                if line.strip()[0] == '|' or line.strip()[-1] == '|':
                    num_lines += 1
                    continue
            
            is_table = False
        
        # since tables take more than one line, increment accordingly
        self.increment = num_lines

        # the number of rows are one line less, due to the extra comment line in markdown
        rows = num_lines - 1

        # get markdown table and convert to HTML
        markdown_table = self.content[index + 1:index + num_lines]
        markdown_converter = Markdown(extras=["tables"])
        html_table = markdown_converter.convert('\n'.join(markdown_table))

        # check if last row has a newline tag, if so, remove it
        if html_table[-1] == '\n':
            html_table = html_table[:-1]

        # add table caption support
        caption_comment = self.caption_extractor.extract(self.content[index], ComponentType.TABLE)

        # add weird WordPress hack where each table is an HTML figure element ... don't ask me, I think Epstein is behind it!
        html_table = '<figure class="wp-block-table">\n' + caption_comment + html_table + '</figure>\n'

        # add table styling with beautiful soup
        soup = BeautifulSoup(html_table, "html.parser")

        for table in soup.find_all("table"):
            table["class"] = "has-fixed-layout"

        for th in soup.find_all("th"):
            th["class"] = "has-text-align-center"
            th["data-align"] = "center"

        for td in soup.find_all("td"):
            td["class"] = "has-text-align-center"
            td["data-align"] = "center"

        return soup.decode(formatter=None)
