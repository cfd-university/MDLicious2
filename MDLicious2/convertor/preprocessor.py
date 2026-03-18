from MDLicious2.components.componentManager import ComponentManager

class Preprocessor:
    def __init__(self, content, components):
        self.content = content
        self.components = components
        self.processed_content = ''

    def preprocess_content(self):
        index = 0
        while index < len(self.content):
            line = self.content[index]

            component_match = False
            for component in self.components:
                if component.match(index):
                    self.processed_content += component.convert(index) + '\n'
                    index = component.next(index)
                    component_match = True
                    break

            if not component_match:
                self.processed_content += line + '\n'
                index += 1
