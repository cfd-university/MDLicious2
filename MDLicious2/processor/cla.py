from argparse import ArgumentParser

class CommandLineArguments:
    def __init__(self):
        self.parser = ArgumentParser(
            prog='MDLicious2',
            description='Convert markdown to HTML using cfd.university specific extensions',
            epilog='Tom-Robin Teschner'
        )

        self.parser.add_argument('-i', '--input', type=str, required=True, help='Input file')
        self.parser.add_argument('-o', '--output', type=str, required=True, help='Output file')

        self.input = self.__get_input()
        self.output = self.__get_output()
    
    def __get_input(self):
        return self.parser.parse_args().input

    def __get_output(self):
        return self.parser.parse_args().output