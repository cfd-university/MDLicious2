import re
from re import sub

from markdown2 import Markdown
from bs4 import BeautifulSoup, NavigableString
from emoji import emojize

from MDLicious2.javascriptRuntime import convert_latex_equation

SHORTCODE_PATTERN = re.compile(
    r'^[ \t]*\[custom_category_posts_list category_slug="[^"]*"\][ \t]*$',
    re.MULTILINE
)

# emoji shortcodes must not be substituted inside code listings, inline code
# spans, or rendered katex markup (katex emits both MathML and HTML)
EMOJI_PROTECTED_TAGS = ('pre', 'code', 'math', 'script', 'style')

class Mark2HTML:
    def __init__(self, content):
        self.content = content
        self.converter = Markdown(extras=["toc"])
        self.protected_shortcodes = []

    def convert(self):
        # needs to run first so markdown does not destroy formating
        self.content = self.__convert_inline_equations(self.content)
        self.content = self.__protect_shortcodes(self.content)

        # convert markdown to html
        html = self.converter.convert(self.content)

        # modify html with custom behaviour
        html = self.__insert_toc(html)
        html = self.__add_security_to_links(html)
        html = self.__add_class_to_blockquotes(html)
        html = self.__convert_emojis(html)
        html = self.__prettify(html)
        html = self.__restore_shortcodes(html)

        return html.split('\n')

    def __insert_toc(self, html):
        temp_html = html.split('\n')
        for i in range(len(temp_html)):
            if temp_html[i].lower().find('[toc]') != -1:
                temp_html[i] = html.toc_html
        html = '\n'.join(temp_html)
        return html

    def __add_security_to_links(self, html):
        soup = BeautifulSoup(html, "html.parser")
        for a in soup.find_all("a", href=True):
            href = a["href"]

            # ignore links to current page
            if href.startswith("#"):
                continue

            # ignore links to cfd.university pages
            if not href.startswith(("http://", "https://")):
                continue
            
            # all other pages get get additional security measures
            # to avoid malicious javascript injection (tabnabbing)
            a["target"] = "_blank"
            a["rel"] = "noopener noreferrer"

        return soup.decode(formatter="minimal")

    def __add_class_to_blockquotes(self, html):
        soup = BeautifulSoup(html, "html.parser")
        for blockquote in soup.find_all("blockquote"):
            blockquote["class"] = "wp-block-quote"

        return soup.decode(formatter="minimal")

    def __convert_emojis(self, html):
        soup = BeautifulSoup(html, "html.parser")

        # only plain text nodes are eligible, comments and other preformatted
        # strings would lose their markers when replaced by a NavigableString
        for text in soup.find_all(string=lambda s: type(s) is NavigableString):
            if self.__is_emoji_protected(text):
                continue

            converted = emojize(str(text), language='alias')
            if converted != text:
                text.replace_with(converted)

        return soup.decode(formatter="minimal")

    def __is_emoji_protected(self, text):
        for parent in text.parents:
            if parent.name in EMOJI_PROTECTED_TAGS:
                return True

            classes = parent.get("class") or []
            if any(css_class.startswith("katex") for css_class in classes):
                return True

        return False

    def __protect_shortcodes(self, content):
        self.protected_shortcodes = []

        def stash(match):
            self.protected_shortcodes.append(match.group(0).strip())
            return f'SHORTCODEPLACEHOLDER{len(self.protected_shortcodes) - 1}'

        return SHORTCODE_PATTERN.sub(stash, content)

    def __restore_shortcodes(self, html):
        for i, shortcode in enumerate(self.protected_shortcodes):
            placeholder = f'SHORTCODEPLACEHOLDER{i}'
            html = html.replace(f'<p>{placeholder}</p>', shortcode)
            html = html.replace(placeholder, shortcode)
        return html

    def __convert_inline_equations(self, html):
        return sub(r'(?<!\$)\$(?!\$)(.*?)(?<!\$)\$(?!\$)', self.replace_inline_katex, html)
    
    def replace_inline_katex(self, match):
        equation = match.group(1)
        return convert_latex_equation(equation, 'false')

    def __prettify(self, html):
        soup = BeautifulSoup(html, "html.parser")

        # bs4 removes new lines, adding them back after block level elements
        elements = ["p", "ul", "ol", "h1", "h2", "h3", "h4", "h5", "h6",
                    "iframe", "img", "div", "figure", "table", "blockquote", "figcaption"]
        for tag in soup.find_all(elements):
            tag.insert_after("\n")

        return soup.decode(formatter="minimal")