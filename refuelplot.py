import seaborn as sns
import matplotlib.style
import matplotlib as mpl
import matplotlib.pyplot as plt


# reFUEL color palette for up to 10 colors
COLORS = (
    '#c72321', '#0d8085', '#efc220', '#fbd7a8', '#1b494d',
    '#861719', '#ba9f7c', '#7a6952', '#6e9b9e', '#af8f19',
)


# TODO this should not be used with German keys I guess, translate first!
TOPIC_COLORS = {
    'Kernkraft':    ('#460e0e', '#7a1718', '#a51f22', '#c42528', '#f8494b'),
    'Erdöl':        ('#3b3329', '#796953', '#ba9f7f', '#dfbf99', '#f9d6ac'),
    'Steinkohle':   ('#20262d', '#323c47', '#445261', '#57697c', '#96a4b4'),
    'Braunkohle':   ('#2d2322', '#463534', '#725f5d', '#957f7d', '#b49b9a'),
    'Wasser':       ('#0f4241', '#136663', '#1e8e88', '#af8f19', '#91cec7'),
    'Wind':         ('#273738', '#246b71', '#6a9395', '#84bcbf', '#9bdade'),
    'Photovoltaik': ('#705b10', '#b08f19', '#d6ae1e', '#f0c322', '#ffde65'),
    'Erdgas':       ('#493614', '#74551a', '#987019', '#d29918', '#e9bb60'),
    'Biomasse':     ('#383c19', '#585e23', '#757e2e', '#757e2e', '#bcc46b'),
}


def setup():
    """Call this function before creating plots to enable reFUEL plotting styles.
    Interesting resources to extend this function include:
    - https://matplotlib.org/tutorials/introductory/customizing.html
    - https://hfstevance.com/blog/2019/7/22/matplotlib-style
    - https://seaborn.pydata.org/tutorial/aesthetics.html
    """
    # TODO future possible parameter of this function: context (as in seaborn)

    sns.set_palette(sns.color_palette(COLORS))
