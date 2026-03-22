from MDLicious2.components.base import Component
from MDLicious2.components.captionExtractor import CaptionExtractor

class Figure(Component):
    def __init__(self, content):
        super().__init__(content)

    def match(self, index):
        line = self.content[index]
        is_figure = line.find('<!-- figure:') != -1
        has_source = False

        if is_figure:
            next_line = self.content[index + 1]
            has_source = next_line.find('![') != -1

        return is_figure and has_source

    def convert(self, index):
        # since figures take two lines, we need to increment the index by two to avoid converting the figure twice
        self.increment = 2

        # get the next two lines for processing
        line = self.content[index]
        next_line = self.content[index + 1]

        # extract meta information from image
        width = line.split('width: ')[1].split('px')[0]
        alt_text = next_line.split('![')[1].split('](')[0]
        image_source = next_line.split('](')[1].split(')')[0]

        caption = CaptionExtractor()

        figure_html = '<figure class="wp-block-image aligncenter size-large is-resized">\n'
        figure_html += f'<img src="{image_source}" alt="{alt_text}" class="wp-image-5550" style="width:{width}px"/>'
        figure_html += caption.extract(line)
        figure_html += '</figure>'

        return figure_html
