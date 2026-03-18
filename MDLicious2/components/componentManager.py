

class ComponentManager:
    def __init__(self):
        self.components = []

    def register(self, component):
        self.components.append(component)