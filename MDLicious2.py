from MDLicious2 import CommandLineArguments, FileProcessor, Mark2HTML, Preprocessor
from MDLicious2 import ComponentManager, YouTube, EquationEnvironment, Code, Figure, Table


def main():
    # process input file
    args = CommandLineArguments()
    filereader = FileProcessor(args)

    # register customised markdown to html converters
    component_manager = ComponentManager()
    component_manager.register(YouTube(filereader.content))
    component_manager.register(EquationEnvironment(filereader.content))
    component_manager.register(Code(filereader.content))
    component_manager.register(Figure(filereader.content))
    component_manager.register(Table(filereader.content))

    # preprocess markdown file with customised converters
    # they will already convert markdown to html
    preprocessor = Preprocessor(filereader.content, component_manager.components)
    preprocessor.preprocess_content()
        
    # convert remaining (standard) markdown to html
    converter = Mark2HTML(preprocessor.processed_content)
    html_content = converter.convert()    

    # write output file
    filereader.output(html_content)


if __name__ == '__main__':
    main()