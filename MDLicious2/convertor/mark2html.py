from markdown2 import Markdown

class Mark2HTML:
    def __init__(self, content):
        self.content = content
        self.converter = Markdown(extras=["toc", "fenced-code-blocks", "tables", "header-ids"])

    def convert(self):
        html = '<html>\n<head>\n<link rel="stylesheet" href="code.css">\n</head>\n<body>\n'
        html += self.converter.convert(self.content)
        html += '\n</body>\n</html>'
        toc = ''

        # html = self.converter.convert(self.content)
        # toc = html.toc_html
        return toc + html
        