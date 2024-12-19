import matplotlib.pyplot as plt
from ..config import install_dir
import os 

def set_style(style='default'):
    style_path = os.path.join(install_dir, 'plotting', 'stylelib', f'{style}.mplstyle')
    assert os.path.exists(style_path), f"Style file {style} not found."
    plt.style.use(style_path)