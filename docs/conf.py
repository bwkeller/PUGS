project = "PUGS"
copyright = "2024, BW Keller"
author = "BW Keller"

extensions = [
    "myst_parser",
    "sphinx_design",
]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

myst_enable_extensions = ["colon_fence", "dollarmath"]

html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "github_url": "https://github.com/bwkeller/PUGS",
    "use_edit_page_button": True,
}
html_context = {
    "github_user": "bwkeller",
    "github_repo": "PUGS",
    "github_version": "main",
    "doc_path": "docs",
}
