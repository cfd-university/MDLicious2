from MDLicious2.components.base import Component

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
        has_tag = False

        for i in range(start, end):
            if self.content[i].find(r'\tag{eq:') != -1:
                has_tag = True
                tag_value = self.content[i].replace(r'\tag{eq:', '').replace('}', '')
                tag = r'\tag{eq:' + tag_value + '}\n'
            else:
                equation += '\n' + self.content[i]
        
        # manage automatic equation tag increase (automatic equation numbering)
        if has_tag:
            equation += tag
        else:
            equation += r'\tag{eq:tag-' + str(self.tag_counter) + '}\n'
        self.tag_counter += 1

        return f'<div class="wp-block-katex-display-block katex-eq" data-katex-display="true"><pre>{equation}</pre></div>'