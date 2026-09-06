project = "PUGS"
copyright = "2024, BW Keller"
author = "BW Keller"

extensions = [
    "myst_parser",
    "sphinx_design",
]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

myst_enable_extensions = ["colon_fence", "dollarmath"]

html_theme = "nasa1976"
html_short_title = "PUGS"
html_theme_options = {
    "program_label": "Portable Universal Galaxy Sampler",
    "document_number": "PUGS",
    "color_scheme": "light",
}
