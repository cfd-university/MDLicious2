from os.path import join
import re

from MDLicious2.components.base import ComponentType

class CaptionMatcher:
    def __init__(self):
        self.counter = {
            ComponentType.FIGURE: 1,
            ComponentType.TABLE: 1,
            ComponentType.CODE: 1,
            ComponentType.EQUATION: 1
        }

        self.counter_map = {}
    
    def setup_equation_tags(self, content):
        index = 0
        equation_counter = 0
        while (index < len(content)):
            line = content[index].strip()
            is_equation = line.strip().find('$$') != -1 and len(line.strip()) == 2

            if is_equation:
                num_lines_equations = 0
                equation_counter += 1
                while True:
                    line = content[index + 1 + num_lines_equations].strip()
                    is_end = line.strip().find('$$') != -1 and len(line.strip()) == 2
                    num_lines_equations += 1

                    if is_end:
                        break
                
                has_tag = False
                tag = ''
                for i in range(index + 1, index + num_lines_equations):
                    line = content[i].strip()
                    if line.find(r'\tag') != -1:
                        has_tag = True
                        tag = self.__extract_tag(line)
                        break
                
                if has_tag:
                    self.counter_map[tag] = self.counter[ComponentType.EQUATION]
                elif not has_tag:
                    tag = f'eq:equation-{equation_counter}'
                    last_eq_line = index + num_lines_equations
                    content.insert(last_eq_line, r'\tag{' + tag + '}')
                    self.counter_map[tag] = self.counter[ComponentType.EQUATION]

                # increment equation counter           
                self.counter[ComponentType.EQUATION] += 1
                index += num_lines_equations + 2
            index += 1

        return content

    def setup_ref_map(self, content):
        for index in range(0, len(content)):
            previous_line = content[index - 1].strip()
            line = content[index].strip()

            if self.__is_figure(line):
                tag = self.__extract_tag(line)
                
                if len(tag) > 0:
                    self.counter_map[tag] = self.counter[ComponentType.FIGURE]
                else:
                    temp = f'fig:figure-{self.counter[ComponentType.FIGURE]}'
                    self.counter_map[temp] = self.counter[ComponentType.FIGURE]
                self.counter[ComponentType.FIGURE] += 1

            elif self.__is_table(line):
                tag = self.__extract_tag(line)
                
                if len(tag) > 0:
                    self.counter_map[tag] = self.counter[ComponentType.TABLE]
                else:
                    temp = f'tab:table-{self.counter[ComponentType.TABLE]}'
                    self.counter_map[temp] = self.counter[ComponentType.TABLE]
                self.counter[ComponentType.TABLE] += 1

            elif self.__is_code(line):
                tag = self.__extract_tag(previous_line)
                
                if len(tag) > 0:
                    self.counter_map[tag] = self.counter[ComponentType.CODE]
                    self.counter[ComponentType.CODE] += 1

            # elif self.__is_equation(line):
            #     tag = self.__extract_tag(line)
                
            #     if len(tag) > 0:
            #         self.counter_map[tag] = self.counter[ComponentType.EQUATION]
            #     else:
            #         temp = f'eq:equation-{self.counter[ComponentType.EQUATION]}'
            #         self.counter_map[temp] = self.counter[ComponentType.EQUATION]
            #     self.counter[ComponentType.EQUATION] += 1

    def __is_figure(self, line):
        return line.find('<!-- figure') != -1

    def __is_table(self, line):
        return line.find('<!-- table') != -1

    def __is_code(self, line):
        return re.search(r'^```[A-Za-z0-9+]+\s*$', line)

    def __is_equation(self, line):
        return line.find(r'\tag{eq:') != -1

    def __has_tag(self, line):
        return line.find(r'\tag') != -1

    def __has_ref(self, line):
        return line.find(r'\ref') != -1

    def __extract_tag(self, line):
        if self.__has_tag(line):
            return line.split(r'\tag')[1].split('{')[1].split('}')[0]
        else:
            return ''
    
    def __extract_ref(self, line):
        if self.__has_ref(line):
            return line.split(r'\ref')[1].split('{')[1].split('}')[0]
        else:
            return ''

    def substitute(self, content):
        content = content.split('\n')
        for index in range(0, len(content)):
            line = content[index].strip()

            if self.__has_ref(line):
                counter = 0
                has_ref = content[index].find(r'\ref') != -1
                while has_ref:
                    for tag in self.counter_map.keys():
                        ref_found = self.__extract_ref(content[index])
                        if ref_found.strip() == tag.strip():
                            content[index] = content[index].replace(r'\ref{' + tag + '}', f'{self.counter_map[tag]}')

                        # if we have a typo in our ref tag, it will never be substituted
                        # if this is the case, silently fail
                        # this will be caught later by the ref checker
                        if counter > 100:
                            has_ref = False
                        counter += 1
            
            # replace equation tag explicitly by number
            if line.find(r'\tag{eq:') != -1:
                for tag in self.counter_map.keys():
                    if content[index].find(tag) != -1:
                        content[index] = content[index].replace(f'{tag}', f'{ self.counter_map[tag]}')

        return '\n'.join(content)
