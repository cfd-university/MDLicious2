import json
from sys import argv

class CommandLineArguments:
    def __init__(self):
        assert len(argv) == 2, "Invalid number of arguments, requires 1 json file as argument!"

        with open(argv[1], 'r') as f:
            self.parser = json.load(f)

        self.input = self.parser["inputFile"]
        self.output = self.parser["outputFile"]
        self.replace = self.parser["replace"]