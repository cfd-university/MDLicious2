

class FileProcessor:
    def __init__(self, cla):
        self.input_file = cla.input
        self.output_file = cla.output

        self.content = self.__read_file()

    def __read_file(self):
        with open(self.input_file) as f:
            return f.read().splitlines()

    def output(self, content):
        with open(self.output_file, 'w') as f:
            f.write(content)