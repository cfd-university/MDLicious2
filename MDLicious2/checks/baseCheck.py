from abc import abstractmethod, ABCMeta

class BaseCheck(metaclass=ABCMeta):
    def __init__(self):
        pass

    @abstractmethod
    def run(self, line, warnings, errors):
        pass

    def _create_warning(self, warning, line):
        return f'<p>{warning}</p><p>{line}</p>'

    def _create_error(self, error, line):
        return f'<p>{error}</p><p>{line}</p>'
    
    def _highlight(self, tag, line):
        return line.replace(tag, f'<mark>{tag}</mark>')
