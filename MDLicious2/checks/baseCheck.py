from abc import abstractmethod, ABCMeta

class BaseCheck(metaclass=ABCMeta):
    def __init__(self):
        pass

    @abstractmethod
    def run(self, line, warnings, errors):
        pass
