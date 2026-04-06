from os.path import join
import re

from MDLicious2.components.base import ComponentType

class CaptionMatcher:
    def __init__(self, content):
        self.content = content.split('\n')

        self.counter = {
            ComponentType.FIGURE: 1,
            ComponentType.TABLE: 1,
            ComponentType.CODE: 1,
            ComponentType.EQUATION: 1
        }

        self.counter_map = {}

        self.__scan()

    def __scan(self):
        for index in range(0, len(self.content)):
            previous_line = self.content[index - 1].strip()
            line = self.content[index].strip()

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

            elif self.__is_equation(line):
                tag = self.__extract_tag(line)
                
                if len(tag) > 0:
                    self.counter_map[tag] = self.counter[ComponentType.EQUATION]
                else:
                    temp = f'eq:equation-{self.counter[ComponentType.EQUATION]}'
                    self.counter_map[temp] = self.counter[ComponentType.EQUATION]
                self.counter[ComponentType.EQUATION] += 1

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

    def substitute(self):
        for index in range(0, len(self.content)):
            line = self.content[index].strip()

            if self.__has_ref(line):
                counter = 0
                has_ref = self.content[index].find(r'\ref') != -1
                while has_ref:
                    for tag in self.counter_map.keys():
                        ref_found = self.__extract_ref(self.content[index])
                        if ref_found.strip() == tag.strip():
                            self.content[index] = self.content[index].replace(r'\ref{' + tag + '}', f'{self.counter_map[tag]}')

                        # if we have a typo in our ref tag, it will never be substituted
                        # if this is the case, silently fail
                        # this will be caught later by the ref checker
                        if counter > 100:
                            has_ref = False
                        counter += 1
            
            # replace equation tag explicitly by number
            if line.find(r'\tag{eq:') != -1:
                for tag in self.counter_map.keys():
                    if self.content[index].find(tag) != -1:
                        self.content[index] = self.content[index].replace(f'{tag}', f'{ self.counter_map[tag]}')

        return '\n'.join(self.content)
