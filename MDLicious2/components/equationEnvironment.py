from MDLicious2.components.base import Component
from MDLicious2.javascriptRuntime import convert_latex_equation

class EquationEnvironment(Component):
    def __init__(self, content):
        super().__init__(content)
        self.tag_counter = 1

    def match(self, index):
        line = self.content[index]
        return (line.find('$$') != -1) and (len(line.strip()) == 2)

    def convert(self, index):
        start, end = self._find_start_end_based_on_pattern(index, '$$')

        equation = ''
        for i in range(start, end):
            equation += '\n' + self.content[i]
        
        return convert_latex_equation(equation, 'true')