from os.path import basename, join, splitext

class FileProcessor:
    def __init__(self, cla):
        self.input_file = cla.input
        self.output_directory = cla.output

        file_name = splitext(basename(self.input_file))[0] + '.html'

        self.output_file = join(self.output_directory, file_name)
        self.replace = cla.replace

        self.content = self.__read_file()
        self.__replace()

    def __read_file(self):
        with open(self.input_file, encoding='utf-8') as f:
            return f.read().splitlines()

    def __replace(self):
        for i in range(0, len(self.content)):
            for key, value in self.replace.items():
                self.content[i] = self.content[i].replace(key, value)

    def output(self, content):
        with open(self.output_file, 'w', encoding='utf-8') as f:
            f.write('\n'.join(content))