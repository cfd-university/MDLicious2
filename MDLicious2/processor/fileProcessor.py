

class FileProcessor:
    def __init__(self, cla):
        self.input_file = cla.input
        self.output_file = cla.output
        self.replace = cla.replace

        self.content = self.__read_file()
        self.__replace()

    def __read_file(self):
        with open(self.input_file) as f:
            return f.read().splitlines()

    def __replace(self):
        for i in range(0, len(self.content)):
            for key, value in self.replace.items():
                self.content[i] = self.content[i].replace(key, value)

    def output(self, content):
        with open(self.output_file, 'w') as f:
            f.write(content)