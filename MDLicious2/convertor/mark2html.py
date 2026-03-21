from markdown2 import Markdown

class Mark2HTML:
    def __init__(self, content):
        self.content = content
        self.converter = Markdown(extras=["toc", "fenced-code-blocks", "tables", "header-ids"])

    def convert(self):
        html = self.converter.convert(self.content)
        html = self.__insert_toc(html)
        return html

    def __insert_toc(self, html):
        temp_html = html.split('\n')
        for i in range(len(temp_html)):
            if temp_html[i].find('[toc]') != -1:
                temp_html[i] = html.toc_html
        html = '\n'.join(temp_html)
        return html   