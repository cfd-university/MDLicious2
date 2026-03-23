from os.path import join
import re

from MDLicious2.components.base import Component
from MDLicious2.convertor.captionExtractor import ComponentType

from pygments import highlight
from pygments.lexers import get_lexer_by_name
from pygments.formatters import HtmlFormatter

class Code(Component):
    def __init__(self, content):
        super().__init__(content)

    def match(self, index):
        line = self.content[index]
        return re.search(r'^```[A-Za-z0-9+]+\s*$', line) 

    def convert(self, index):
        start, end = self._find_start_end_based_on_pattern(index, '```')

        # extract possible caption on the previous line
        caption = self.caption_extractor.extract(self.content[index - 1], ComponentType.CODE)

        language = 'text'
        if len(self.content[index].strip()) > 3:
            language = self.content[index].strip()[3:]

        lexer = get_lexer_by_name(language)
        formatter = HtmlFormatter(style='nord', linenos="table", cssclass="codehilite")
        css = formatter.get_style_defs('.codehilite')

        # # uncomment the following to generate new code styles
        # with open('code.css', 'w') as f:
        #     f.write(css)

        code = ''
        for i in range(start, end):
            code += '\n' + self.content[i]
        code += '\n'

        html_code = highlight(code, lexer, formatter)

        return html_code + caption