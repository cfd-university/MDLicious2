from re import sub

from MDLicious2.components.base import Component

class InlineEquation(Component):
    def __init__(self, content):
        super().__init__(content)

    def match(self, index):
        line = self.content[index]
        has_inline_equation = (line.find('$') != -1) and (line.find('$$') == -1)
        return has_inline_equation

    def convert(self, index):
        line = self.content[index]
        return sub(r'(?<!\$)\$(?!\$)(.*?)(?<!\$)\$(?!\$)', r'[katex]\1[/katex]', line)