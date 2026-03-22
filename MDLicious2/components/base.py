from abc import abstractmethod, ABCMeta

class Component(metaclass=ABCMeta):
    def __init__(self, content):
        self.content = content
        self.increment = 1

    @abstractmethod
    def match(self, index):
        pass

    @abstractmethod
    def convert(self, index):
        pass

    def next(self, index):
        return index + self.increment

    def _find_start_end_based_on_pattern(self, index, pattern):
        self.__set_increment_based_on_pattern(index, pattern)
        start_index = self.__multiline_start_index(index)
        end_index   = self.__multiline_end_index(index)

        return start_index, end_index

    def _set_increment(self, increment):
        self.increment = increment

    def __set_increment_based_on_pattern(self, index, end_pattern):
        start_index = index + 1
        for temp_index in range(start_index, len(self.content)):
            if self.content[temp_index].find(end_pattern) != -1:
                self.increment = temp_index - start_index + 2
                return
    
    def __multiline_start_index(self, index):
        return index + 1
    
    def __multiline_end_index(self, index):
        return index + self.increment - 1