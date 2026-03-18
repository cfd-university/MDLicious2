from MDLicious2.components.base import Component

class YouTube(Component):
    def __init__(self, content):
        super().__init__(content)

    def match(self, index):
        line = self.content[index]
        return line.find('<!-- youtube: https://www.youtube.com/watch?v=') != -1


    def convert(self, index):
        line = self.content[index]
        url = line.replace('<!-- youtube: https://www.youtube.com/watch?v=', '')
        url = url.replace(' -->', '')

        if url.find('&') != -1:
            url = url.split('&')[0]
        
        return f'<iframe src="https://www.youtube-nocookie.com/embed/{url}" title="YouTube video player" frameborder="0" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" referrerpolicy="strict-origin-when-cross-origin" allowfullscreen style="width:100%; aspect-ratio:16/9; border:0;"></iframe>'